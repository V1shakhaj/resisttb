"""Train per-drug resistance models on CRyPTIC TB genomes and compare them with the catalogue.

Inputs (from the CRyPTIC June 2022 release, see README.md):
  EFFECTS.csv.gz           mutations each sample has in the resistance genes for each drug,
                           plus the catalogue's R/S/U call for that mutation
  UKMYC_PHENOTYPES.csv.gz  lab-measured resistance (R/S) for each sample and drug

Outputs (in --out):
  metrics.json  held-out test performance of the model and the catalogue, per drug
  model.json    learned mutation weights per drug, used by the web demo
  report.md     plain-language summary of the results
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

DRUGS = {
    "RIF": "Rifampicin", "INH": "Isoniazid", "EMB": "Ethambutol", "LEV": "Levofloxacin",
    "MXF": "Moxifloxacin", "AMI": "Amikacin", "KAN": "Kanamycin", "ETH": "Ethionamide",
    "RFB": "Rifabutin", "BDQ": "Bedaquiline", "LZD": "Linezolid", "DLM": "Delamanid", "CFZ": "Clofazimine",
}
SEED = 42
LOW_DATA = 200  # fewer resistant samples than this is flagged as unreliable in the report
MIN_SAMPLES = 3  # a mutation must appear in at least this many training samples to become a feature


def load_effects(path):
    """Read the default-catalogue rows only, in chunks so it fits on a small cloud machine."""
    cols = ["UNIQUEID", "DRUG", "GENE", "MUTATION", "PREDICTION", "DEFAULT_CATALOGUE"]
    parts = []
    for chunk in pd.read_csv(path, usecols=cols, chunksize=2_000_000, low_memory=False):
        chunk = chunk[(chunk.DEFAULT_CATALOGUE == True) & chunk.DRUG.isin(DRUGS)]  # noqa: E712
        parts.append(chunk.drop(columns="DEFAULT_CATALOGUE"))
    df = pd.concat(parts, ignore_index=True)
    # GARC grammar: a trailing x/z/o marks a null, mixed or minor-population call. Keep confident calls only.
    df = df[~df.MUTATION.str[-1].str.lower().isin(["x", "z", "o"])]
    df["FEATURE"] = df.GENE + "@" + df.MUTATION
    return df


def load_phenotypes(path):
    p = pd.read_csv(path, low_memory=False)
    p = p[p.DRUG.isin(DRUGS) & p.BINARY_PHENOTYPE.isin(["R", "S"]) & p.PHENOTYPE_QUALITY.isin(["HIGH", "MEDIUM"])]
    return p[["UNIQUEID", "DRUG", "BINARY_PHENOTYPE"]]


def sens_spec(y_true, y_pred):
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    tp = int(((y_true == 1) & (y_pred == 1)).sum()); fn = int(((y_true == 1) & (y_pred == 0)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum()); fp = int(((y_true == 0) & (y_pred == 1)).sum())
    return {
        "sensitivity": tp / (tp + fn) if tp + fn else None,
        "specificity": tn / (tn + fp) if tn + fp else None,
        "tp": tp, "fn": fn, "tn": tn, "fp": fp,
    }


def train_drug(drug, effects, pheno):
    labels = pheno[pheno.DRUG == drug].set_index("UNIQUEID").BINARY_PHENOTYPE.eq("R").astype(int)
    eff = effects[(effects.DRUG == drug) & effects.UNIQUEID.isin(labels.index)]
    # Only samples that were sequenced (appear somewhere in EFFECTS) can be scored
    labels = labels[labels.index.isin(effects.UNIQUEID.unique())]
    n_r = int(labels.sum())
    if n_r < 20:
        return None, {"drug": drug, "name": DRUGS[drug], "skipped": f"only {n_r} resistant samples"}

    ids = labels.index.to_numpy()
    train_ids, test_ids = train_test_split(ids, test_size=0.2, random_state=SEED, stratify=labels.values)
    train_set = set(train_ids)

    counts = eff[eff.UNIQUEID.isin(train_set)].groupby("FEATURE").UNIQUEID.nunique()
    features = sorted(counts[counts >= MIN_SAMPLES].index)
    fidx = {f: i for i, f in enumerate(features)}
    row = {u: i for i, u in enumerate(ids)}
    e = eff[eff.FEATURE.isin(fidx)].drop_duplicates(["UNIQUEID", "FEATURE"])
    X = csr_matrix((np.ones(len(e)), (e.UNIQUEID.map(row), e.FEATURE.map(fidx))), shape=(len(ids), len(features)))
    y = labels.loc[ids].to_numpy()
    tr = np.array([row[u] for u in train_ids]); te = np.array([row[u] for u in test_ids])

    model = LogisticRegression(penalty="l1", solver="liblinear", C=1.0, class_weight="balanced", max_iter=2000)
    model.fit(X[tr], y[tr])
    prob = model.predict_proba(X[te])[:, 1]
    ml = sens_spec(y[te], (prob >= 0.5).astype(int))
    ml["auc"] = float(roc_auc_score(y[te], prob))

    # Catalogue baseline: resistant if any mutation is graded R; U (unknown) counts as susceptible
    calls = eff[eff.UNIQUEID.isin(set(test_ids))].groupby("UNIQUEID").PREDICTION.agg(set)
    cat_pred = np.array([1 if "R" in calls.get(u, set()) else 0 for u in test_ids])
    has_u = np.array(["U" in calls.get(u, set()) and "R" not in calls.get(u, set()) for u in test_ids])
    cat = sens_spec(y[te], cat_pred)
    cat["unknown_only_fraction"] = float(has_u.mean())

    weights = {f: round(float(w), 3) for f, w in zip(features, model.coef_[0]) if abs(w) > 1e-6}
    top = sorted(weights.items(), key=lambda kv: -kv[1])[:10]
    metrics = {
        "drug": drug, "name": DRUGS[drug], "n_samples": len(ids), "n_resistant": n_r,
        "n_test": len(te), "n_features": len(features), "n_features_used": len(weights),
        "model": ml, "catalogue": cat, "top_resistance_mutations": top,
    }
    return {"intercept": round(float(model.intercept_[0]), 3), "weights": weights}, metrics


def pct(x):
    return "n/a" if x is None else f"{100 * x:.1f}%"


def write_report(metrics, path):
    lines = [
        "# TB drug-resistance model: results", "",
        "Trained on the CRyPTIC dataset (lab-tested TB isolates with whole-genome sequences).",
        "Each drug has its own model. 80% of samples were used for training and 20% were held back for testing.",
        "Numbers below are on the held-back 20% only.", "",
        "- **Sensitivity**: of the truly resistant samples, how many the method caught. Missing these is the dangerous error.",
        "- **Specificity**: of the truly susceptible samples, how many were correctly called susceptible.",
        "- **Catalogue**: the rule-based approach (known resistance mutations), with unknown mutations treated as susceptible.", "",
        "| Drug | Samples (resistant) | Model sensitivity | Model specificity | Catalogue sensitivity | Catalogue specificity |",
        "|---|---|---|---|---|---|",
    ]
    for m in metrics:
        if "skipped" in m:
            lines.append(f"| {m['name']} | skipped: {m['skipped']} | | | | |")
            continue
        cat = m["catalogue"]
        no_rules = cat["tp"] + cat["fp"] == 0
        cat_cols = "no rules | no rules" if no_rules else f"{pct(cat['sensitivity'])} | {pct(cat['specificity'])}"
        flag = " ⚠️" if m["n_resistant"] < LOW_DATA else ""
        lines.append(
            f"| {m['name']}{flag} | {m['n_samples']:,} ({m['n_resistant']:,}) | {pct(m['model']['sensitivity'])} | "
            f"{pct(m['model']['specificity'])} | {cat_cols} |"
        )
    lines += [
        "",
        f"⚠️ Fewer than {LOW_DATA} resistant samples to learn from. Treat these results as unreliable until more data is added.",
        "",
        "\"No rules\" means the catalogue used in this dataset makes no resistance calls for that drug, so the model is the only predictor.",
        "",
        "## Strongest resistance signals the model learned",
        "",
    ]
    for m in metrics:
        if "skipped" not in m:
            muts = ", ".join(f"`{f.replace('@', ' ')}`" for f, _ in m["top_resistance_mutations"][:5])
            lines.append(f"- **{m['name']}**: {muts}")
    lines += ["", "Research prototype. Not validated for clinical use."]
    Path(path).write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    print("Loading phenotypes..."); pheno = load_phenotypes(Path(args.data) / "UKMYC_PHENOTYPES.csv.gz")
    print("Loading mutations (takes a few minutes)..."); effects = load_effects(Path(args.data) / "EFFECTS.csv.gz")

    models, metrics = {}, []
    for drug in DRUGS:
        model, m = train_drug(drug, effects, pheno)
        metrics.append(m)
        if model:
            models[drug] = model
            print(f"{drug}: model sens {pct(m['model']['sensitivity'])} spec {pct(m['model']['specificity'])} | "
                  f"catalogue sens {pct(m['catalogue']['sensitivity'])} spec {pct(m['catalogue']['specificity'])}")
        else:
            print(f"{drug}: {m['skipped']}")

    (out / "metrics.json").write_text(json.dumps(metrics, indent=2))
    (out / "model.json").write_text(json.dumps({"drugs": DRUGS, "models": models}, separators=(",", ":")))
    write_report(metrics, out / "report.md")
    print(f"Wrote {out}/metrics.json, model.json, report.md")


if __name__ == "__main__":
    main()
