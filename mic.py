"""Predict how resistant a TB sample is: its MIC (minimum inhibitory concentration), not just yes/no.

A yes/no tool such as TB-Profiler says "resistant" or "susceptible". For dosing it also matters how far
above the cut-off a sample sits: low-level resistance may still respond to a higher dose.

Model: one ridge regression per drug on the doubling-dilution scale. Each mutation shifts the predicted
MIC up or down by a learned number of doublings, so every prediction can be explained.

Validation showed the model clearly beats a yes/no call on resistant samples only for isoniazid and
kanamycin. Kanamycin is no longer recommended by WHO, so the website predicts a level only for isoniazid. For every drug it also shows the typical MIC range
measured in samples carrying each resistance mutation (a lookup, not a prediction).

Accuracy is reported as essential agreement (predicted MIC within one doubling of the lab's), the
standard measure for MIC methods, and compared with two baselines on the same samples:
  - most common MIC: always predict the drug's most common MIC
  - WHO yes/no call: knowing only the WHO catalogue's resistant/susceptible call (what TB-Profiler
    reports), predict the typical MIC of training samples with that call

Inputs (in --data): UKMYC_PHENOTYPES.csv.gz, EFFECTS.csv.gz, SITES.csv, WHO_catalogue_v2.1.csv
Also reads results/model.json (from train.py) to know which mutations are linked to resistance.
Outputs (in --out): mic_model.json, mic_metrics.json, mic_report.md
"""
import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split

from evaluate import WhoCatalogue, country_of
from train import DRUGS, MIN_SAMPLES, SEED, load_effects, pct

LEVEL_DRUGS = ["INH"]  # drugs where the level model clearly beat the yes/no baseline and is still in use
LOW_LEVEL_STEPS = 2  # up to this many doublings above the cut-off counts as low-level resistance
STRONG_WEIGHT = 1.0  # resistance-model weight above which a mutation counts as resistance-linked
MIN_MUTATION_SAMPLES = 10  # solo samples needed before showing a typical MIC range for a mutation


def parse_mic(s):
    s = str(s)
    if s.startswith("<="):
        return float(s[2:]), "le"
    if s.startswith(">"):
        return float(s[1:]), "gt"
    return float(s), "eq"


def mic_scale(rows):
    """Put both plate designs on one doubling-dilution scale and return each sample's bin.

    The plates test slightly different ranges, so the scale runs from the highest bottom concentration
    to the lowest top concentration that both plates share. Values outside it fall into the end bins.
    """
    parsed = rows.MIC.map(parse_mic)
    rows = rows.assign(CONC=parsed.str[0], KIND=parsed.str[1])
    lo = rows[rows.KIND == "le"].groupby("PLATEDESIGN").CONC.min().max()
    hi = rows[rows.KIND == "gt"].groupby("PLATEDESIGN").CONC.max().min()
    steps = int(round(np.log2(hi / lo)))
    ladder = [lo * 2 ** k for k in range(steps + 1)]  # lo, 2*lo, ..., hi
    labels = [f"≤{lo:g}"] + [f"{c:g}" for c in sorted(rows[(rows.KIND == "eq") & (rows.CONC > lo * 1.01) & (rows.CONC <= hi * 1.01)].CONC.unique())] + [f">{hi:g}"]
    if len(labels) != steps + 2:  # fall back to a computed ladder if a dilution is never observed
        labels = [f"≤{lo:g}"] + [f"{c:.3g}" for c in ladder[1:]] + [f">{hi:g}"]

    def to_bin(conc, kind):
        if kind == "le" or conc <= lo * 1.01:
            return 0
        if kind == "gt" or conc > hi * 1.01:
            return steps + 1
        return int(round(np.log2(conc / lo)))

    bins = np.array([to_bin(c, k) for c, k in zip(rows.CONC, rows.KIND)])
    return bins, labels


def matrix(fsets, ids, idx):
    rows, cols = [], []
    for r, u in enumerate(ids):
        for f in fsets.get(u, ()):
            if f in idx:
                rows.append(r); cols.append(idx[f])
    return csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(ids), len(idx)))


def fit(fsets, ids, y):
    counts = Counter(f for u in ids for f in fsets.get(u, ()))
    idx = {f: i for i, f in enumerate(sorted(f for f, c in counts.items() if c >= MIN_SAMPLES))}
    model = Ridge(alpha=1.0).fit(matrix(fsets, ids, idx), y)
    return model, idx


def predict(model, idx, fsets, ids, top):
    return np.clip(np.rint(model.predict(matrix(fsets, ids, idx))), 0, top).astype(int)


def agreement(true, pred, ecoff):
    true, pred = np.asarray(true), np.asarray(pred)
    tr, pr = true > ecoff, pred > ecoff
    return {
        "essential_agreement": float(np.mean(np.abs(true - pred) <= 1)),
        # the same, for truly resistant samples only: the ones where the level matters
        "essential_agreement_resistant": float(np.mean(np.abs(true[tr] - pred[tr]) <= 1)) if tr.any() else None,
        "exact_agreement": float(np.mean(true == pred)),
        "category_agreement": float(np.mean(tr == pr)),
        # missed resistance (dangerous) and false resistance, as a share of truly resistant / susceptible samples
        "very_major_error": float(np.mean(~pr[tr])) if tr.any() else None,
        "major_error": float(np.mean(pr[~tr])) if (~tr).any() else None,
    }


def evaluate_split(fsets, who_r, ids, y, train, test, ecoff, top):
    model, idx = fit(fsets, ids[train], y[train])
    out = {"model": agreement(y[test], predict(model, idx, fsets, ids[test], top), ecoff)}
    mode = Counter(y[train]).most_common(1)[0][0]
    out["most_common_mic"] = agreement(y[test], np.full(test.sum(), mode), ecoff)
    if who_r is not None:
        # Typical MIC for samples the WHO catalogue calls resistant / not resistant
        med = {flag: int(np.median(y[train][who_r[train] == flag])) if (who_r[train] == flag).any() else mode for flag in (True, False)}
        out["who_yes_no"] = agreement(y[test], np.array([med[f] for f in who_r[test]]), ecoff)
    else:
        out["who_yes_no"] = None
    return out


def mutation_mic(fsets, ids, y, labels, positive, ecoff):
    """Typical MIC (median and middle half) for each resistance-linked mutation.

    Uses only samples where it is the sole resistance-linked mutation for the drug ("solo", as in the WHO
    catalogue method), so mutations that often travel together don't borrow each other's MIC.
    Only mutations whose typical MIC is above the cut-off are kept; the rest are strain markers.
    """
    by = {}
    for u, v in zip(ids, y):
        linked = [f for f in fsets.get(u, ()) if f in positive]
        if len(linked) == 1:
            by.setdefault(linked[0], []).append(v)
    out = {}
    for f, v in by.items():
        if len(v) >= MIN_MUTATION_SAMPLES:
            q1, med, q3 = np.percentile(v, [25, 50, 75], method="nearest")
            if med <= ecoff:
                continue
            out[f] = {"n": len(v), "median": labels[int(med)], "low": labels[int(q1)], "high": labels[int(q3)]}
    return out


def run_drug(drug, pheno, effects, who, sites, positive):
    rows = pheno[pheno.DRUG == drug]
    y, labels = mic_scale(rows)
    ids = rows.UNIQUEID.to_numpy()
    top = len(labels) - 1
    # ECOFF: the highest MIC bin where most samples are labelled susceptible
    s_frac = pd.Series(rows.BINARY_PHENOTYPE.eq("S").to_numpy()).groupby(y).mean()
    ecoff = int(max(b for b, f in s_frac.items() if f >= 0.5))

    eff = effects[effects.DRUG == drug]
    fsets = eff.groupby("UNIQUEID").FEATURE.agg(set).to_dict()
    all_feats = effects[effects.UNIQUEID.isin(set(ids))].groupby("UNIQUEID").FEATURE.agg(set).to_dict()
    who_r = (np.array(["R" in who.sample_calls(all_feats.get(u, ()), drug) for u in ids]) if drug in who.drugs else None)

    res = {"drug": drug, "name": DRUGS[drug], "n": len(ids), "bins": labels, "ecoff_bin": ecoff, "ecoff": labels[ecoff]}
    tr_ids, te_ids = train_test_split(np.arange(len(ids)), test_size=0.2, random_state=SEED, stratify=y > ecoff)
    train = np.zeros(len(ids), bool); train[tr_ids] = True
    res["random_split"] = evaluate_split(fsets, who_r, ids, y, train, ~train, ecoff, top)

    countries = np.array([country_of(u, sites) for u in ids])
    preds, true, who_pred, mode_pred = [], [], [], []
    for c in np.unique(countries):
        test = countries == c
        model, idx = fit(fsets, ids[~test], y[~test])
        preds.append(predict(model, idx, fsets, ids[test], top)); true.append(y[test])
        mode = Counter(y[~test]).most_common(1)[0][0]
        mode_pred.append(np.full(test.sum(), mode))
        if who_r is not None:
            med = {f: int(np.median(y[~test][who_r[~test] == f])) if (who_r[~test] == f).any() else mode for f in (True, False)}
            who_pred.append(np.array([med[f] for f in who_r[test]]))
    true = np.concatenate(true)
    res["unseen_countries"] = {
        "model": agreement(true, np.concatenate(preds), ecoff),
        "most_common_mic": agreement(true, np.concatenate(mode_pred), ecoff),
        "who_yes_no": agreement(true, np.concatenate(who_pred), ecoff) if who_r is not None else None,
    }

    model_json = {"bins": labels, "ecoff_bin": ecoff, "mutations": mutation_mic(fsets, ids, y, labels, positive, ecoff)}
    if drug in LEVEL_DRUGS:
        final, idx = fit(fsets, ids, y)
        model_json["intercept"] = round(float(final.intercept_), 3)
        model_json["weights"] = {f: round(float(w), 3) for f, w in zip(sorted(idx, key=idx.get), final.coef_) if abs(w) >= 0.05}
    res["level_on_website"] = drug in LEVEL_DRUGS
    return res, model_json


def write_report(results, path):
    def row(r, test):
        t = r[test]
        who, m = t["who_yes_no"], t["model"]
        name = f"**{r['name']}** ✓" if r["level_on_website"] else r["name"]
        return (f"| {name} | {r['n']:,} | {r['ecoff']} | {pct(m['essential_agreement_resistant'])} | "
                f"{pct(who['essential_agreement_resistant']) if who else 'no rules'} | {pct(m['essential_agreement'])} | "
                f"{pct(who['essential_agreement']) if who else 'no rules'} | {pct(t['most_common_mic']['essential_agreement'])} |")
    head = ["| Drug | Samples | Cut-off (ECOFF, mg/L) | Resistant samples: AI model | Resistant samples: WHO yes/no | All samples: AI model | All samples: WHO yes/no | All samples: most common MIC |",
            "|---|---|---|---|---|---|---|---|"]
    lines = [
        "# MIC prediction: how resistant, not just whether", "",
        "The MIC is the lowest drug concentration that stops the bacteria growing. Samples just above the cut-off",
        "(low-level resistance) may still respond to a higher dose; samples far above it will not.", "",
        "**Essential agreement**: the share of samples where the predicted MIC is within one doubling of the lab's.",
        "This is the standard accuracy measure for MIC methods. Higher is better.", "",
        "- **AI model**: ridge regression on mutations, one model per drug.",
        "- **WHO yes/no call**: knowing only the WHO catalogue's resistant/susceptible call (what TB-Profiler reports),",
        "  the best guess is the typical MIC for that call. The gap between this column and the AI model is the value",
        "  of predicting the level, not just yes or no.",
        "- **Most common MIC**: always guessing the drug's most common MIC. Any useful method must beat this.",
        "- **Resistant samples** columns count only samples above the cut-off, where the level matters for treatment.", "",
        "## Result", "",
        "On resistant samples, the AI model clearly beats the yes/no call only for **isoniazid**, where it separates high-level",
        "*katG* from low-level *inhA* resistance, and for kanamycin, which WHO no longer recommends. It roughly ties for most other",
        "drugs and is worse for moxifloxacin. Isoniazid (✓) is therefore the only drug where the website predicts a level.",
        "For other drugs, knowing which resistance mutation is present already says most of what there is to know about the",
        "level, so the website shows the typical MIC range for each mutation instead. It is measured in samples where that is",
        f"the only resistance-linked mutation for the drug (at least {MIN_MUTATION_SAMPLES} such samples).",
        f"Low-level resistance means up to {LOW_LEVEL_STEPS} doublings above the cut-off; high-level means further above.", "",
        "## Unseen countries (each country predicted by a model trained on the others)", "", *head,
        *[row(r, "unseen_countries") for r in results], "",
        "## Random 80/20 split", "", *head, *[row(r, "random_split") for r in results], "",
        "## Caveats", "",
        "- MICs were measured on CRyPTIC UKMYC plates. Other lab methods give somewhat different values.",
        "- The two plate versions test different ranges, so values beyond the shared range are grouped in the end bins (≤ and >).",
        "- Effects of mutations are treated as additive on the doubling scale, which is a simplification.",
        "- Dosing decisions need a specialist and local guidance. Research prototype, not for clinical use.",
    ]
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    data, out = Path(args.data), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    pheno = pd.read_csv(data / "UKMYC_PHENOTYPES.csv.gz", low_memory=False)
    pheno = pheno[pheno.DRUG.isin(DRUGS) & pheno.PHENOTYPE_QUALITY.isin(["HIGH", "MEDIUM"]) & pheno.MIC.notna()]
    effects = load_effects(data / "EFFECTS.csv.gz")
    pheno = pheno[pheno.UNIQUEID.isin(set(effects.UNIQUEID))]  # sequenced samples only
    sites_df = pd.read_csv(data / "SITES.csv", dtype=str)
    sites = dict(zip(sites_df.SITEID, sites_df.COUNTRY.fillna("Other")))
    who = WhoCatalogue(data / "WHO_catalogue_v2.1.csv")

    clf = json.loads((out / "model.json").read_text())["models"]
    results, models = [], {}
    for drug in DRUGS:
        positive = {f for f, w in clf.get(drug, {}).get("weights", {}).items() if w >= STRONG_WEIGHT}
        res, model = run_drug(drug, pheno, effects, who, sites, positive)
        results.append(res); models[drug] = model
        u = res["unseen_countries"]
        print(f"{drug}: resistant samples within one doubling (unseen countries): model "
              f"{pct(u['model']['essential_agreement_resistant'])}, WHO yes/no "
              f"{pct(u['who_yes_no']['essential_agreement_resistant']) if u['who_yes_no'] else 'n/a'}; "
              f"{len(model['mutations'])} mutations with typical MIC")

    (out / "mic_metrics.json").write_text(json.dumps(results, indent=2))
    (out / "mic_model.json").write_text(json.dumps({"low_level_steps": LOW_LEVEL_STEPS, "drugs": models}, separators=(",", ":")))
    write_report(results, out / "mic_report.md")
    print(f"Wrote {out}/mic_metrics.json, mic_model.json, mic_report.md")


if __name__ == "__main__":
    main()
