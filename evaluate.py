"""Validate the models on countries and studies they were not trained on, and benchmark them against
the WHO mutation catalogue (2023, v2.1) and the CRyPTIC catalogue on exactly the same samples.

Two tests:
  1. Unseen countries: for each country in the CRyPTIC plate data, train on all other countries and
     predict that country. Every sample is predicted by a model that never saw its country.
  2. Other studies: train on all CRyPTIC plate data, then predict samples from 7 other studies
     (MGIT, LJ and other lab methods, 40+ countries). Samples from patients in the training data are removed.

The WHO catalogue is applied with piezo, the reference rule engine for this mutation grammar.
TB-Profiler bases its resistance calls on the same WHO catalogue, so this is the catalogue-level
equivalent of a TB-Profiler benchmark.

Extra inputs (in --data, see README.md): DST_MEASUREMENTS.csv.gz, SITES.csv, WHO_catalogue_v2.1.csv
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

from train import DRUGS, LOW_DATA, MIN_SAMPLES, load_effects, load_phenotypes, new_model, pct, sens_spec

MAIN_DRUGS = ["RIF", "INH", "EMB", "LEV", "MXF", "AMI"]  # shown per country / per study in the report
MIN_GROUP = 20  # per-country and per-study rows need at least this many resistant and susceptible samples


def subject(uid):
    return ".".join(uid.split(".")[:4])  # site.XX.subj.YYY: one patient


def country_of(uid, sites):
    return sites.get(uid.split(".")[1], "Other")


def load_external(path, exclude_subjects):
    d = pd.read_csv(path, dtype=str)
    d = d[(d.SOURCE != "CRyPTIC") & d.DRUG.isin(DRUGS) & d.PHENOTYPE.isin(["R", "S"])]
    d = d[~d.UNIQUEID.map(subject).isin(exclude_subjects)]
    # One label per sample and drug. Drop pairs where the lab results disagree.
    g = d.groupby(["UNIQUEID", "DRUG"]).agg(N=("PHENOTYPE", "nunique"), PHENOTYPE=("PHENOTYPE", "first"), SOURCE=("SOURCE", "first"))
    return g[g.N == 1].drop(columns="N").reset_index()


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


def evaluate_drug(drug, effects, pheno, external, sites, who, all_feats):
    eff = effects[effects.DRUG == drug]
    fsets = eff.groupby("UNIQUEID").FEATURE.agg(set).to_dict()
    cryptic_sets = eff.groupby("UNIQUEID").PREDICTION.agg(set).to_dict()
    cryptic_has_rules = (eff.PREDICTION == "R").any()
    cryptic = (lambda u: cryptic_sets.get(u, set())) if cryptic_has_rules else None
    who_fn = (lambda u: who.sample_calls(all_feats.get(u, ()), drug)) if drug in who.drugs else None

    sequenced = set(all_feats)
    lab = pheno[(pheno.DRUG == drug) & pheno.UNIQUEID.isin(sequenced)]
    ids = lab.UNIQUEID.to_numpy()
    y = lab.BINARY_PHENOTYPE.eq("R").astype(int).to_numpy()
    result = {"drug": drug, "name": DRUGS[drug], "low_data": int(y.sum()) < LOW_DATA}

    # 1. Leave one country out
    countries = np.array([country_of(u, sites) for u in ids])
    prob = np.zeros(len(ids))
    for c in np.unique(countries):
        test = countries == c
        if y[~test].sum() == 0:
            continue
        prob[test] = fit_predict(fsets, ids[~test], y[~test], ids[test])
    result["unseen_countries"] = score(ids, y, prob, cryptic, who_fn)
    per_country = {}
    for c in np.unique(countries):
        m = countries == c
        if MIN_GROUP <= y[m].sum() <= m.sum() - MIN_GROUP:
            per_country[c] = score(ids[m], y[m], prob[m], cryptic, who_fn)
    result["by_country"] = per_country

    # 2. Other studies: train on everything from CRyPTIC plates
    ext = external[(external.DRUG == drug) & external.UNIQUEID.isin(sequenced)]
    if len(ext) and ext.PHENOTYPE.eq("R").any():
        eids = ext.UNIQUEID.to_numpy()
        ey = ext.PHENOTYPE.eq("R").astype(int).to_numpy()
        eprob = fit_predict(fsets, ids, y, eids)
        result["other_studies"] = score(eids, ey, eprob, cryptic, who_fn)
        by_study = {}
        for s in sorted(ext.SOURCE.unique()):
            m = ext.SOURCE.to_numpy() == s
            if MIN_GROUP <= ey[m].sum() <= m.sum() - MIN_GROUP:
                by_study[s] = score(eids[m], ey[m], eprob[m], cryptic, who_fn)
        result["by_study"] = by_study
    else:
        result["other_studies"], result["by_study"] = None, {}
    return result


def cell(r, key):
    if r is None or r.get(key) is None:
        return "no rules | no rules"
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
        flag = " ⚠️" if r["low_data"] else ""
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
        "## Test 1: countries the model never saw", "",
        "For each country, the model was trained on all other countries and then tested on that one.", "",
        *comparison_table(results, "unseen_countries"), "",
        f"### By country (at least {MIN_GROUP} resistant and {MIN_GROUP} susceptible samples)", "",
        *group_table(results, "by_country", "Country"), "",
        "## Test 2: other studies", "",
        "The model was trained on all CRyPTIC plate data and tested on samples from 7 other studies. These use different lab methods",
        "(mostly MGIT liquid culture and solid media at WHO critical concentrations) and come from over 40 countries.",
        "Patients who also appear in the CRyPTIC plate data were removed.", "",
        *comparison_table(results, "other_studies"), "",
        f"### By study (at least {MIN_GROUP} resistant and {MIN_GROUP} susceptible samples)", "",
        *group_table(results, "by_study", "Study"), "",
        "## Caveats", "",
        "- **The WHO catalogue had a head start.** WHO built its 2023 catalogue from about 52,000 isolates, including the CRyPTIC data",
        "  and several of the studies tested here. The catalogue has partly seen these samples before, so its numbers here are optimistic.",
        "  The AI model never saw any test sample, so its numbers are the fairer estimate for new data.",
        f"- ⚠️ Fewer than {LOW_DATA} resistant samples in the training data. Treat as unreliable.",
        "- The mutation table records insertions and deletions without their exact length, so WHO rules for frameshifts cannot fire.",
        "  This mainly lowers WHO catalogue sensitivity for bedaquiline, clofazimine and delamanid.",
        "- Only mutations in genes the CRyPTIC catalogue links to a drug are available. WHO rules for other genes are missed.",
        "- A true TB-Profiler run needs raw sequencing reads, several terabytes for these samples, which this project does not process.",
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
    args = ap.parse_args()
    data, out = Path(args.data), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    pheno = load_phenotypes(data / "UKMYC_PHENOTYPES.csv.gz")
    effects = load_effects(data / "EFFECTS.csv.gz")
    sites_df = pd.read_csv(data / "SITES.csv", dtype=str)
    sites = dict(zip(sites_df.SITEID, sites_df.COUNTRY.fillna("Other")))
    external = load_external(data / "DST_MEASUREMENTS.csv.gz", set(pheno.UNIQUEID.map(subject)))
    who = WhoCatalogue(data / "WHO_catalogue_v2.1.csv")

    needed = set(pheno.UNIQUEID) | set(external.UNIQUEID)
    all_feats = effects[effects.UNIQUEID.isin(needed)].groupby("UNIQUEID").FEATURE.agg(set).to_dict()
    print(f"{len(all_feats):,} sequenced samples with lab results "
          f"({external.UNIQUEID.isin(all_feats.keys()).sum():,} lab results from other studies)")

    results = []
    for drug in DRUGS:
        r = evaluate_drug(drug, effects, pheno, external, sites, who, all_feats)
        results.append(r)
        u, o = r["unseen_countries"], r["other_studies"]
        line = f"{drug}: unseen countries model sens {pct(u['model']['sensitivity'])} spec {pct(u['model']['specificity'])}"
        if o:
            w = o["who_catalogue"]
            line += (f" | other studies n={o['n']:,} model sens {pct(o['model']['sensitivity'])} spec {pct(o['model']['specificity'])}"
                     + (f", WHO sens {pct(w['sensitivity'])} spec {pct(w['specificity'])}" if w else ", WHO no rules"))
        print(line)

    (out / "validation.json").write_text(json.dumps(results, indent=2))
    write_report(results, out / "validation.md")
    print(f"Wrote {out}/validation.json, validation.md")


if __name__ == "__main__":
    main()
