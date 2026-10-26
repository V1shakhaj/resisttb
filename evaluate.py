"""Validate the models on countries and studies they were not trained on, and benchmark them against
the WHO mutation catalogue (2023, v2.1) and the CRyPTIC catalogue on exactly the same samples.

The model is trained on the CRyPTIC plate data plus 7 other studies. Each test trains it on everything except
the group being tested:
  1. Unseen countries: each country's CRyPTIC plate samples, predicted by a model trained on all other countries
     (from every study).
  2. Unseen studies: each of the 7 other studies (MGIT, LJ and other lab methods, 40+ countries), predicted by a
     model trained on CRyPTIC and the other 6 studies.
  3. Unseen lineages: each main strain lineage's CRyPTIC plate samples, predicted by a model trained on all other
     lineages. Lineages differ in background mutations, so this checks the model isn't relying on strain markers.
Patients in the end-to-end TB-Profiler check are never used.

The WHO catalogue is applied with piezo, the reference rule engine for this mutation grammar.
TB-Profiler bases its resistance calls on the same WHO catalogue, so this is the catalogue-level
equivalent of a TB-Profiler benchmark.

Extra inputs (in --data, see README.md): DST_MEASUREMENTS.csv.gz, DST_SAMPLES.csv.gz, SITES.csv, MYKROBE_LINEAGE.csv.gz,
WHO_catalogue_v2.1.csv
Outputs (in --out): validation.json, validation.md
"""
import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import piezo
from scipy.sparse import csr_matrix
from sklearn.metrics import roc_auc_score

from train import (DRUGS, LOW_DATA, MIN_SAMPLES, held_out_ids, load_effects, load_external, model_filters, model_rows,
                   new_model, pct, sens_spec, subject, training_labels)

MAIN_DRUGS = ["RIF", "INH", "EMB", "LEV", "MXF", "AMI"]  # shown per country / per study in the report
MAIN_LINEAGES = ["Lineage 1", "Lineage 2", "Lineage 3", "Lineage 4"]
MIN_GROUP = 20  # per-country and per-study rows need at least this many resistant and susceptible samples


def country_of(uid, sites):
    return sites.get(uid.split(".")[1], "Other")


class WhoCatalogue:
    def __init__(self, path):
        self.cat = piezo.ResistanceCatalogue(str(path))
        self.drugs = set(self.cat.catalogue.drugs)
        self.cache = {}

    def calls(self, feature):
        if feature not in self.cache:
            try:
                r = self.cat.predict(feature)
            except Exception:
                r = None
            self.cache[feature] = r if isinstance(r, dict) else {}
        return self.cache[feature]

    def sample_calls(self, features, drug):
        return {self.calls(f).get(drug, "S") for f in features}


def fit_predict(fsets, train_ids, y_train, test_ids):
    counts = Counter(f for u in train_ids for f in fsets.get(u, ()))
    idx = {f: i for i, f in enumerate(sorted(f for f, c in counts.items() if c >= MIN_SAMPLES))}

    def matrix(ids):
        rows, cols = [], []
        for r, u in enumerate(ids):
            for f in fsets.get(u, ()):
                if f in idx:
                    rows.append(r); cols.append(idx[f])
        return csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(ids), len(idx)))

    return new_model().fit(matrix(train_ids), y_train).predict_proba(matrix(test_ids))[:, 1]


def score(ids, y, prob, cryptic, who):
    out = {"n": len(ids), "n_resistant": int(y.sum())}
    if prob is None:  # no model for drugs with too little data
        out["model"] = None
    else:
        out["model"] = sens_spec(y, (prob >= 0.5).astype(int))
        out["model"]["auc"] = float(roc_auc_score(y, prob)) if 0 < y.sum() < len(y) else None
    for name, calls in (("who_catalogue", who), ("cryptic_catalogue", cryptic)):
        if calls is None:
            out[name] = None
            continue
        sets = [calls(u) for u in ids]
        r = sens_spec(y, np.array([1 if "R" in s else 0 for s in sets]))
        r["uncertain_fraction"] = float(np.mean([("U" in s and "R" not in s) for s in sets])) if sets else None
        out[name] = r
    return out


def evaluate_drug(drug, effects, pool, countries, lineages, who, all_feats, filters):
    eff = effects[effects.DRUG == drug]
    fsets = model_rows(effects, drug, filters).groupby("UNIQUEID").FEATURE.agg(set).to_dict()
    cryptic_sets = eff.groupby("UNIQUEID").PREDICTION.agg(set).to_dict()
    cryptic_has_rules = (eff.PREDICTION == "R").any()
    cryptic = (lambda u: cryptic_sets.get(u, set())) if cryptic_has_rules else None
    who_fn = (lambda u: who.sample_calls(all_feats.get(u, ()), drug)) if drug in who.drugs else None

    lab = pool[(pool.DRUG == drug) & pool.UNIQUEID.isin(set(all_feats))]
    ids, y, src = lab.UNIQUEID.to_numpy(), lab.R.to_numpy(), lab.SOURCE.to_numpy()
    plate = src == "CRyPTIC"
    low = int(y.sum()) < LOW_DATA
    result = {"drug": drug, "name": DRUGS[drug], "low_data": low}
    sub = lambda arr, m: None if arr is None else arr[m]

    def held_out_test(groups, test_mask, keys):
        """Predict each group in test_mask with a model trained on every sample outside that group."""
        if low:
            return None
        prob = np.full(len(ids), np.nan)
        for k in keys:
            test = test_mask & (groups == k)
            train = groups != k
            if test.any() and y[train].sum():
                prob[test] = fit_predict(fsets, ids[train], y[train], ids[test])
        return prob

    def tested(prob, mask):
        return mask if prob is None else mask & ~np.isnan(prob)

    # 1. Unseen countries: CRyPTIC plate samples, by country
    ctry = np.array([countries.get(u, "Unknown") for u in ids])
    cprob = held_out_test(ctry, plate, np.unique(ctry[plate]))
    m = tested(cprob, plate)
    result["unseen_countries"] = score(ids[m], y[m], sub(cprob, m), cryptic, who_fn)
    result["by_country"] = {}
    for c in np.unique(ctry[m]):
        g = m & (ctry == c)
        if MIN_GROUP <= y[g].sum() <= g.sum() - MIN_GROUP:
            result["by_country"][c] = score(ids[g], y[g], sub(cprob, g), cryptic, who_fn)

    # 3. Unseen lineages: CRyPTIC plate samples of the four main lineages
    lin = np.array([lineages.get(u, "Unknown") for u in ids])
    lprob = held_out_test(lin, plate & np.isin(lin, MAIN_LINEAGES), MAIN_LINEAGES)
    m = tested(lprob, plate & np.isin(lin, MAIN_LINEAGES))
    result["unseen_lineages"] = score(ids[m], y[m], sub(lprob, m), cryptic, who_fn)
    result["by_lineage"] = {}
    for l in MAIN_LINEAGES:
        g = m & (lin == l)
        if MIN_GROUP <= y[g].sum() <= g.sum() - MIN_GROUP:
            result["by_lineage"][l] = score(ids[g], y[g], sub(lprob, g), cryptic, who_fn)

    # 2. Unseen studies: each of the other studies, predicted by a model trained on CRyPTIC and the rest
    studies = [st for st in np.unique(src) if st != "CRyPTIC"]
    sprob = held_out_test(src, ~plate, studies)
    m = tested(sprob, ~plate)
    if m.any() and y[m].any():
        result["other_studies"] = score(ids[m], y[m], sub(sprob, m), cryptic, who_fn)
        result["by_study"] = {}
        for st in studies:
            g = m & (src == st)
            if MIN_GROUP <= y[g].sum() <= g.sum() - MIN_GROUP:
                result["by_study"][st] = score(ids[g], y[g], sub(sprob, g), cryptic, who_fn)
    else:
        result["other_studies"], result["by_study"] = None, {}
    return result


def cell(r, key):
    if r is None or r.get(key) is None:
        return "too little data | too little data" if key == "model" else "no rules | no rules"
    return f"{pct(r[key]['sensitivity'])} | {pct(r[key]['specificity'])}"


def comparison_table(results, test):
    lines = [
        "| Drug | Samples (resistant) | AI model sens. | AI model spec. | WHO catalogue sens. | WHO catalogue spec. | CRyPTIC catalogue sens. | CRyPTIC catalogue spec. |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        t = r.get(test)
        if not t:
            lines.append(f"| {r['name']} | no lab results | | | | | | |")
            continue
        flag = ""
        lines.append(f"| {r['name']}{flag} | {t['n']:,} ({t['n_resistant']:,}) | {cell(t, 'model')} | {cell(t, 'who_catalogue')} | {cell(t, 'cryptic_catalogue')} |")
    return lines


def group_table(results, key, label):
    lines = [f"| {label} | Drug | Samples (resistant) | AI model sens. / spec. | WHO catalogue sens. / spec. |", "|---|---|---|---|---|"]
    for r in results:
        if r["drug"] not in MAIN_DRUGS:
            continue
        for g, t in r[key].items():
            who = t["who_catalogue"]
            w = f"{pct(who['sensitivity'])} / {pct(who['specificity'])}" if who else "no rules"
            lines.append(f"| {g} | {r['name']} | {t['n']:,} ({t['n_resistant']:,}) | "
                         f"{pct(t['model']['sensitivity'])} / {pct(t['model']['specificity'])} | {w} |")
    return lines


def write_report(results, path):
    lines = [
        "# Validation: unseen countries, other studies, and the WHO catalogue", "",
        "These tests check whether the AI model still works on data it was not built from,",
        "and compare it with the rule-based catalogues on exactly the same samples.", "",
        "- **Sensitivity (sens.)**: of the truly resistant samples, how many were caught. Missing these is the dangerous error.",
        "- **Specificity (spec.)**: of the truly susceptible samples, how many were correctly called susceptible.",
        "- **WHO catalogue**: WHO mutation catalogue, 2nd edition (2023, v2.1), applied with piezo. "
        "TB-Profiler bases its resistance calls on this catalogue, so this column is the catalogue-level stand-in for TB-Profiler.",
        "- **CRyPTIC catalogue**: the catalogue shipped with the CRyPTIC data (v1.31).",
        "- Catalogue mutations graded *uncertain* are counted as susceptible, as a clinical report would.", "",
        "The model is trained on the CRyPTIC plate data plus 7 other studies (about 39,000 patients). Each test trains it",
        "on everything except the group being tested. Patients in the end-to-end TB-Profiler check are never used.", "",
        "## Test 1: countries the model never saw", "",
        "Each country's CRyPTIC plate samples were predicted by a model trained on all other countries, from every study.", "",
        *comparison_table(results, "unseen_countries"), "",
        f"### By country (at least {MIN_GROUP} resistant and {MIN_GROUP} susceptible samples)", "",
        *group_table(results, "by_country", "Country"), "",
        "## Test 2: studies the model never saw", "",
        "Each of the 7 other studies was predicted by a model trained on CRyPTIC and the other 6 studies. These studies use",
        "different lab methods (mostly MGIT liquid culture and solid media at WHO critical concentrations) and come from over",
        "40 countries. Patients who also appear in the CRyPTIC plate data were removed.", "",
        *comparison_table(results, "other_studies"), "",
        f"### By study (at least {MIN_GROUP} resistant and {MIN_GROUP} susceptible samples)", "",
        *group_table(results, "by_study", "Study"), "",
        "## Test 3: strain lineages the model never saw", "",
        "TB strains fall into lineages (families) that carry different harmless background mutations. Each of the four main",
        "lineages' CRyPTIC plate samples was predicted by a model trained on all other lineages, from every study. If it had",
        "learned strain markers instead of real resistance mutations, it would fail here. Samples of other or unknown lineage",
        "were used for training only.", "",
        *comparison_table(results, "unseen_lineages"), "",
        f"### By lineage (at least {MIN_GROUP} resistant and {MIN_GROUP} susceptible samples)", "",
        "Lineage 2 includes the Beijing family, which is common in East Asia and often drug-resistant. Lineage 4 is the most",
        "widespread worldwide.", "",
        *group_table(results, "by_lineage", "Lineage"), "",
        "## Caveats", "",
        "- The AI model only uses mutations in the genes the WHO catalogue links to each drug, minus mutations the catalogue",
        "  grades as harmless. An earlier version used every candidate gene and failed this lineage test: it had learned strain",
        "  markers and mutations that merely travel with multidrug resistance (amikacin specificity on unseen lineages was 61%,",
        "  kanamycin 35%, and levofloxacin 91% because of the strain marker gyrA G668D).",
        "- **The WHO catalogue had a head start.** WHO built its 2023 catalogue from about 52,000 isolates, including the CRyPTIC data",
        "  and several of the studies tested here. The catalogue has partly seen these samples before, so its numbers here are optimistic.",
        "  The AI model never saw any test sample, so its numbers are the fairer estimate for new data.",
        f"- Drugs with fewer than {LOW_DATA} resistant samples (bedaquiline, linezolid, delamanid, clofazimine) get no AI model.",
        "  Earlier versions trained one anyway and it learned noise; the catalogue columns still show how the rules perform.",
        "- The mutation table records insertions and deletions without their exact length, so WHO rules for frameshifts cannot fire.",
        "  This mainly lowers WHO catalogue sensitivity for bedaquiline, clofazimine and delamanid.",
        "- Only mutations in genes the CRyPTIC catalogue links to a drug are available. WHO rules for other genes are missed.",
        "- For the real TB-Profiler on raw reads, see the end-to-end check (tbprofiler_benchmark.md).",
        "- Different studies use different lab methods and cut-offs, so some disagreement with the lab result is expected for any method.",
        "- When the AI model and the WHO catalogue both disagree with the lab for the same group (for example Sweden),",
        "  the cause is more likely that group's lab results (method, cut-off or sample selection) than either predictor.",
        "", "Research prototype. Not validated for clinical use.",
    ]
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="results")
    ap.add_argument("--held-out", default="results/tbprofiler_samples.json")
    args = ap.parse_args()
    data, out = Path(args.data), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    pool = training_labels(data, held_out_ids(args.held_out))
    effects = load_effects(data / "EFFECTS.csv.gz")
    # Country: from the site for CRyPTIC samples, from the sample record for the other studies (3-letter codes)
    sites_df = pd.read_csv(data / "SITES.csv", dtype=str)
    site_code = dict(zip(sites_df.SITEID, sites_df.COUNTRY_3_LETTER_CODE))
    names = dict(zip(sites_df.COUNTRY_3_LETTER_CODE.dropna(), sites_df.COUNTRY.dropna()))
    samples = pd.read_csv(data / "DST_SAMPLES.csv.gz", dtype=str).drop_duplicates("UNIQUEID").set_index("UNIQUEID")
    countries = {}
    for u in pool.UNIQUEID.unique():
        code = site_code.get(u.split(".")[1]) or (samples.COUNTRY_CODE.get(u) if u in samples.index else None)
        countries[u] = names.get(code, code) if isinstance(code, str) else "Unknown"
    lin_df = pd.read_csv(data / "MYKROBE_LINEAGE.csv.gz", dtype=str)
    lineages = dict(zip(samples.index, samples.MYKROBE_LINEAGE_NAME_1.fillna("Unknown")))
    lineages.update(zip(lin_df.UNIQUEID, lin_df.MYKROBE_LINEAGE_NAME_1.fillna("Unknown")))
    who = WhoCatalogue(data / "WHO_catalogue_v2.1.csv")
    filters = model_filters(data / "WHO_catalogue_v2.1.csv")

    all_feats = effects[effects.UNIQUEID.isin(set(pool.UNIQUEID))].groupby("UNIQUEID").FEATURE.agg(set).to_dict()
    print(f"{len(all_feats):,} sequenced patients with lab results")

    results = []
    for drug in DRUGS:
        r = evaluate_drug(drug, effects, pool, countries, lineages, who, all_feats, filters)
        results.append(r)
        u, o, l = r["unseen_countries"], r["other_studies"], r["unseen_lineages"]
        if u["model"] is None:
            print(f"{drug}: too little data for a model"); continue
        f = lambda t: f"{pct(t['model']['sensitivity'])}/{pct(t['model']['specificity'])}" if t and t["model"] else "n/a"
        print(f"{drug}: countries {f(u)} | lineages {f(l)} | studies {f(o)}")

    (out / "validation.json").write_text(json.dumps(results, indent=2))
    write_report(results, out / "validation.md")
    print(f"Wrote {out}/validation.json, validation.md")


if __name__ == "__main__":
    main()
