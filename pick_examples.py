"""Pick real, anonymised patients for the web demo from the "other studies" test set.

Every patient picked here comes from a study the model was not trained on, and patients who also appear
in the CRyPTIC plate data are excluded, so the model has never seen them.

Patients are chosen by resistance profile (from the lab results) and by how many drugs the lab tested,
never by whether the model gets them right. One extra patient is chosen where the model misses
rifampicin resistance, to show a real failure case.

Inputs (in --data): EFFECTS.csv.gz, UKMYC_PHENOTYPES.csv.gz, DST_MEASUREMENTS.csv.gz, DST_SAMPLES.csv.gz
Uses results/model.json for the miss example.
Output (in --out): examples.json
"""
import argparse
import json
import math
import re
from pathlib import Path

import pandas as pd

from evaluate import subject
from train import load_effects, load_phenotypes

# Demo drug code -> code in the lab results
DEMO_DRUGS = {"RIF": "RIF", "INH": "INH", "EMB": "EMB", "PZA": "PZA", "LFX": "LEV", "MFX": "MXF",
              "BDQ": "BDQ", "LZD": "LZD", "AMK": "AMI", "ETO": "ETH"}
MAX_MUTATIONS = 40  # keep examples readable; very long lists usually mean a mixed or low-quality sample


def display(gene, mutation):
    """GARC spelling -> the gene + change format the demo accepts (it converts back for the model)."""
    m = re.fullmatch(r"([acgt])(-?\d+)([acgt])", mutation)
    if m:
        prefix = "n." if gene in ("rrs", "rrl") else "c."
        return f"{gene} {prefix}{m[2]}{m[1].upper()}>{m[3].upper()}"
    return f"{gene} {mutation.replace('!', '*')}"


def model_prob(model, drug, features):
    md = model["models"].get(drug)
    if not md:
        return None
    z = md["intercept"] + sum(md["weights"].get(f, 0) for f in features)
    return 1 / (1 + math.exp(-z))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    data, out = Path(args.data), Path(args.out)

    pheno = load_phenotypes(data / "UKMYC_PHENOTYPES.csv.gz")
    seen = set(pheno.UNIQUEID.map(subject))
    lab = pd.read_csv(data / "DST_MEASUREMENTS.csv.gz", dtype=str)
    lab = lab[(lab.SOURCE != "CRyPTIC") & lab.PHENOTYPE.isin(["R", "S"]) & ~lab.UNIQUEID.map(subject).isin(seen)]
    wide = lab.groupby(["UNIQUEID", "DRUG"]).PHENOTYPE.agg(lambda x: x.iloc[0] if x.nunique() == 1 else None).unstack()
    info = pd.read_csv(data / "DST_SAMPLES.csv.gz", dtype=str).drop_duplicates("UNIQUEID").set_index("UNIQUEID")
    method = lab.groupby("UNIQUEID").METHOD_2.agg(lambda x: x.mode().iloc[0] if x.notna().any() else "")

    # All mutations in the genes the demo's drugs use
    effects = load_effects(data / "EFFECTS.csv.gz")
    effects = effects[effects.DRUG.isin(set(DEMO_DRUGS.values())) & effects.UNIQUEID.isin(wide.index)]
    muts = effects.drop_duplicates(["UNIQUEID", "FEATURE"]).groupby("UNIQUEID")[["GENE", "MUTATION", "FEATURE"]].agg(list)
    muts = muts[muts.FEATURE.str.len() <= MAX_MUTATIONS]
    wide = wide.loc[wide.index.intersection(muts.index)]

    col = lambda d: wide[d] if d in wide else pd.Series(index=wide.index, dtype=object)
    R, S = (lambda d: col(d).eq("R")), (lambda d: col(d).eq("S"))
    fq_r, fq_s = R("LEV") | R("MXF"), (S("LEV") | S("MXF")) & ~(R("LEV") | R("MXF"))
    grp_a_r = R("BDQ") | R("LZD")
    tested = wide[[c for c in DEMO_DRUGS.values() if c in wide]].notna().sum(axis=1)

    profiles = [
        ("Drug-susceptible", S("RIF") & S("INH") & S("EMB") & fq_s),
        ("MDR-TB", R("RIF") & R("INH") & fq_s),
        ("Pre-XDR", R("RIF") & fq_r & ~grp_a_r),
        ("XDR", R("RIF") & fq_r & grp_a_r),
    ]
    model = json.loads((out / "model.json").read_text())
    rif_prob = pd.Series({u: model_prob(model, "RIF", muts.loc[u, "FEATURE"]) for u in wide.index})
    profiles.append(("Model miss: rifampicin", R("RIF") & S("INH") & (rif_prob < 0.5)))

    examples, used = [], set()
    for label, mask in profiles:
        cands = tested[mask & ~wide.index.isin(used)]
        if cands.empty:
            print(f"{label}: no candidates")
            continue
        # Most drugs tested first; ties broken by ID so the pick is reproducible
        uid = sorted(cands.index, key=lambda u: (-cands[u], u))[0]
        used.add(uid)
        row = info.loc[uid] if uid in info.index else None
        accession = uid.split(".")[3]  # the published, anonymised sample ID (ENA accession for most studies)
        examples.append({
            "id": accession,
            "name": label,
            "study": row["SOURCE"] if row is not None else "",
            "country": row["COUNTRY_CODE"] if row is not None and isinstance(row["COUNTRY_CODE"], str) else "",
            "lab_method": method.get(uid, ""),
            "muts": sorted(display(g, m) for g, m in zip(muts.loc[uid, "GENE"], muts.loc[uid, "MUTATION"])),
            "lab": {demo: wide.loc[uid, code] for demo, code in DEMO_DRUGS.items()
                    if code in wide and isinstance(wide.loc[uid, code], str)},
        })
        print(f"{label}: {accession} ({examples[-1]['study']}, {cands[uid]} drugs tested, {len(examples[-1]['muts'])} mutations)")

    (out / "examples.json").write_text(json.dumps(examples, indent=1))
    print(f"Wrote {out}/examples.json")


if __name__ == "__main__":
    main()
