"""End-to-end check with the real TB-Profiler: raw reads -> TB-Profiler -> ResistTB, compared with lab results.

  select   Pick patients from the "other studies" set (never used for training) that have lab results and public
           paired-end Illumina reads on ENA, mixing resistance profiles. Writes the sample list and job shards.
  compare  Read the TB-Profiler results.json files, apply ResistTB's model to TB-Profiler's own mutation calls,
           and compare TB-Profiler alone, ResistTB's model alone, and TB-Profiler + ResistTB with the lab.

TB-Profiler itself runs in .github/workflows/tbprofiler.yml between these two steps.
"""
import argparse
import json
import math
import re
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from evaluate import load_external, subject
from train import DRUGS, load_phenotypes, pct, sens_spec

SEED = 7
BENCH_DRUGS = ["RIF", "INH", "EMB", "PZA", "LEV", "MXF", "AMI", "KAN", "ETH"]
TBP_DRUG = {"rifampicin": "RIF", "isoniazid": "INH", "ethambutol": "EMB", "pyrazinamide": "PZA", "levofloxacin": "LEV",
            "moxifloxacin": "MXF", "amikacin": "AMI", "kanamycin": "KAN", "ethionamide": "ETH"}
MAX_FASTQ_BYTES = 1_500_000_000  # skip very deep runs to keep each job inside GitHub's limits
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run&format=json&fields=run_accession,fastq_ftp,fastq_bytes,instrument_platform,library_layout"
ACCESSION = re.compile(r"^(SAM[NED][A-Z]?\d+|[ESD]RS\d+)$")


# ---------- select

def ena_runs(accession):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(ENA.format(accession), timeout=30) as r:
                return json.loads(r.read() or b"[]")
        except Exception:
            time.sleep(2 * (attempt + 1))
    return []


def paired_run(accession):
    """The first paired-end Illumina run with both FASTQ files, or None."""
    for run in ena_runs(accession):
        if run.get("instrument_platform") != "ILLUMINA" or run.get("library_layout") != "PAIRED":
            continue
        files = (run.get("fastq_ftp") or "").split(";")
        sizes = [int(b) for b in (run.get("fastq_bytes") or "").split(";") if b.isdigit()]
        r1 = [f for f in files if f.endswith("_1.fastq.gz")]
        r2 = [f for f in files if f.endswith("_2.fastq.gz")]
        if r1 and r2 and sum(sizes) <= MAX_FASTQ_BYTES:
            return {"run": run["run_accession"], "fastq": ["https://" + r1[0], "https://" + r2[0]], "bytes": sum(sizes)}
    return None


def select(args):
    data = Path(args.data)
    pheno = load_phenotypes(data / "UKMYC_PHENOTYPES.csv.gz")
    ext = load_external(data / "DST_MEASUREMENTS.csv.gz", set(pheno.UNIQUEID.map(subject)))
    ext = ext[ext.DRUG.isin(BENCH_DRUGS)]
    wide = ext.pivot(index="UNIQUEID", columns="DRUG", values="PHENOTYPE")
    wide = wide[wide.RIF.notna() & wide.INH.notna()]
    acc = wide.index.str.split(".").str[3]
    wide = wide[(wide.index.str.split(".").str[1] == "ENA") & acc.str.match(ACCESSION.pattern)]
    fq_r = wide.get("LEV", pd.Series(dtype=object)).eq("R") | wide.get("MXF", pd.Series(dtype=object)).eq("R")
    groups = {
        "pre-XDR": wide.RIF.eq("R") & fq_r,
        "RR/MDR": wide.RIF.eq("R") & ~fq_r,
        "Hr-TB": wide.RIF.eq("S") & wide.INH.eq("R"),
        "susceptible": wide.RIF.eq("S") & wide.INH.eq("S"),
    }
    share = {"pre-XDR": 0.15, "RR/MDR": 0.30, "Hr-TB": 0.15, "susceptible": 0.40}
    rng = np.random.default_rng(SEED)
    chosen = []
    for name, mask in groups.items():
        want = round(args.n * share[name])
        pool = list(wide.index[mask])
        rng.shuffle(pool)
        got = 0
        for uid in pool:
            if got >= want:
                break
            run = paired_run(uid.split(".")[3])
            if not run:
                continue
            lab = {d: wide.loc[uid, d] for d in BENCH_DRUGS if d in wide and isinstance(wide.loc[uid, d], str)}
            chosen.append({"id": uid.split(".")[3], "uniqueid": uid, "group": name, "lab": lab, **run})
            got += 1
        print(f"{name}: {got} of {want} (pool {len(pool):,})")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "samples.json").write_text(json.dumps(chosen, indent=1))
    shards = [chosen[i:i + args.per_shard] for i in range(0, len(chosen), args.per_shard)]
    for i, shard in enumerate(shards):
        (out / f"shard_{i}.txt").write_text("".join(f"{s['id']}\t{s['fastq'][0]}\t{s['fastq'][1]}\n" for s in shard))
    (out / "matrix.json").write_text(json.dumps({"shard": list(range(len(shards)))}))
    print(f"Selected {len(chosen)} samples in {len(shards)} shards; {sum(s['bytes'] for s in chosen) / 1e9:.1f} GB of reads")


# ---------- compare

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
       "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
       "Ter": "!", "*": "!"}


def to_feature(v):
    """A TB-Profiler variant in the mutation spelling the model was trained on (same rules as the website)."""
    gene = v.get("gene_name") or v.get("gene") or v.get("locus_tag")
    typ = str(v.get("type") or "").lower()
    if not gene or "synonymous" in typ:
        return None
    ch = v.get("change") or v.get("protein_change") or v.get("nucleotide_change") or ""
    nuc = v.get("nucleotide_change") or ""
    m = re.fullmatch(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)", ch)
    if m and m[1] in AA3 and m[3] in AA3:
        return f"{gene}@{AA3[m[1]]}{m[2]}{AA3[m[3]]}"
    m = re.fullmatch(r"[cn]\.(-?\d+)([ACGT])>([ACGT])", ch)
    if m:
        return f"{gene}@{m[2].lower()}{m[1]}{m[3].lower()}"
    m = re.match(r"c\.(-?\d+)(?:_-?\d+)?(?:del|ins|dup)", nuc)
    if m and re.search(r"fs|frameshift|inframe|ins|del|dup", ch + " " + typ, re.I):
        return f"{gene}@{m[1]}_indel"
    return f"{gene}@{ch}"


def tbp_calls(res):
    calls = {}
    for v in res.get("dr_variants", []):
        for d in v.get("drugs") or v.get("annotation") or []:
            if d.get("type", "drug_resistance") == "drug_resistance":
                code = TBP_DRUG.get(str(d.get("drug", "")).lower())
                if code:
                    calls[code] = True
    return calls


def model_prob(model, drug, features):
    md = model["models"].get(drug)
    if not md:
        return None
    return 1 / (1 + math.exp(-(md["intercept"] + sum(md["weights"].get(f, 0) for f in features))))


def compare(args):
    samples = json.loads(Path(args.samples).read_text())
    model = json.loads(Path(args.model).read_text())
    rows = []
    for s in samples:
        path = Path(args.profiles) / f"{s['id']}.results.json"
        if not path.exists():
            continue
        res = json.loads(path.read_text())
        variants = res.get("dr_variants", []) + res.get("other_variants", [])
        feats = {f for f in map(to_feature, variants) if f}
        calls = tbp_calls(res)
        for drug, lab in s["lab"].items():
            p = model_prob(model, drug, feats)
            rows.append({"id": s["id"], "group": s["group"], "drug": drug, "lab": lab, "tbprofiler": bool(calls.get(drug)),
                         "model_p": p, "lineage": res.get("sub_lineage") or res.get("main_lineage") or ""})
    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    n_done = df.id.nunique() if len(df) else 0
    results = []
    for drug in BENCH_DRUGS:
        d = df[df.drug == drug] if len(df) else df
        if not len(d) or not d.lab.eq("R").any():
            continue
        y = d.lab.eq("R").astype(int).to_numpy()
        tbp = d.tbprofiler.astype(int).to_numpy()
        r = {"drug": drug, "name": DRUGS.get(drug, "Pyrazinamide"), "n": len(d), "n_resistant": int(y.sum()),
             "tbprofiler": sens_spec(y, tbp)}
        if d.model_p.notna().all():
            m = (d.model_p.to_numpy() >= 0.5).astype(int)
            r["resisttb_model"] = sens_spec(y, m)
            r["tbprofiler_plus_resisttb"] = sens_spec(y, np.maximum(tbp, m))
        results.append(r)
    (out / "tbprofiler_benchmark.json").write_text(json.dumps({"samples_profiled": n_done, "samples_selected": len(samples), "drugs": results}, indent=2))
    df.to_csv(out / "tbprofiler_benchmark_calls.csv", index=False)

    def cell(x):
        return f"{pct(x['sensitivity'])} / {pct(x['specificity'])}" if x else "no model"
    lines = [
        "# End-to-end check with the real TB-Profiler", "",
        f"Raw sequencing reads for {n_done} patients (of {len(samples)} selected) were downloaded from ENA and run through",
        "TB-Profiler. ResistTB's model was then applied to TB-Profiler's own mutation calls. None of these patients were",
        "used to train ResistTB. Each cell is sensitivity / specificity against the lab result.", "",
        "- **TB-Profiler**: resistant if TB-Profiler reports a resistance mutation for the drug.",
        "- **ResistTB model**: resistant if ResistTB's model gives 50% or more.",
        "- **TB-Profiler + ResistTB**: resistant if either says so. This is how the add-on works on the website.", "",
        "| Drug | Samples (resistant) | TB-Profiler | ResistTB model | TB-Profiler + ResistTB |",
        "|---|---|---|---|---|",
        *[f"| {r['name']} | {r['n']} ({r['n_resistant']}) | {cell(r['tbprofiler'])} | {cell(r.get('resisttb_model'))} | {cell(r.get('tbprofiler_plus_resisttb'))} |"
          for r in results],
        "", "## Caveats", "",
        "- About 100 patients is a small sample: a difference of one or two patients moves these numbers by several percent.",
        "- Patients were chosen to mix resistance profiles, so the share of resistant samples is higher than in most clinics.",
        "- The lab results come from several studies and methods (mostly MGIT), which have their own error rates.",
        "- TB-Profiler's database and ResistTB's model both draw partly on the WHO catalogue, which was built from some of these studies.",
        "", "Research prototype. Not validated for clinical use.",
    ]
    (out / "tbprofiler_benchmark.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[10:10 + len(results) + 2]))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--data", default="data")
    s.add_argument("--n", type=int, default=100)
    s.add_argument("--per-shard", type=int, default=10)
    s.add_argument("--out", default="benchmark")
    c = sub.add_parser("compare")
    c.add_argument("--samples", default="benchmark/samples.json")
    c.add_argument("--profiles", default="tbp_results")
    c.add_argument("--model", default="results/model.json")
    c.add_argument("--out", default="results")
    args = ap.parse_args()
    select(args) if args.cmd == "select" else compare(args)


if __name__ == "__main__":
    main()
