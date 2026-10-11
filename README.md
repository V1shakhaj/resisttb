# ResistTB Predictor

Predicts drug-resistant tuberculosis from the mutations found when a TB sample is sequenced.
Research prototype. Not for clinical use.

- **Live demo:** https://v1shakhaj.github.io/resisttb/ (`index.html`, served by GitHub Pages)
- **Training code:** `train.py`
- **Validation code:** `evaluate.py`
- **Resistance level (MIC):** `mic.py` predicts how resistant, not just whether ([`results/mic_report.md`](results/mic_report.md))
- **Demo patients:** `pick_examples.py` picks real, anonymised patients the model never saw (`results/examples.json`)
- **Latest results:** [`results/report.md`](results/report.md) (random split) and [`results/validation.md`](results/validation.md) (unseen countries, other studies, WHO catalogue benchmark)

## How it works

1. **Data:** the [CRyPTIC](https://www.crypticproject.org/) June 2022 release: about 12,000 TB isolates with whole-genome sequences and lab-measured resistance to 13 drugs (UKMYC microtitre plates).
   - `EFFECTS.csv.gz`: each sample's mutations in the genes linked to each drug, plus the catalogue's call for each mutation.
   - `UKMYC_PHENOTYPES.csv.gz`: the lab result (resistant or susceptible) per sample and drug. Only high- and medium-quality readings are used.
2. **Model:** one logistic regression per drug with an L1 penalty, which keeps only the mutations that matter. Each mutation gets a weight, so every prediction can be explained. Each drug's model only sees mutations in the genes the WHO catalogue links to that drug; without this, the lineage test showed it learning strain markers.
3. **Evaluation:** a stratified 80/20 split. The 20% test set is never used for training. Results are compared with the catalogue (known resistance mutations, with unknown mutations counted as susceptible).
4. **Validation** (`evaluate.py`):
   - **Unseen countries:** for each country, train on all the others and test on that one.
   - **Other studies:** train on all CRyPTIC plate data, then test on about 40,000 samples from 7 other studies (`DST_MEASUREMENTS.csv.gz`), which use different lab methods (MGIT, LJ) and come from over 40 countries. Patients who appear in the training data are removed.
   - **Unseen lineages:** for each of the four main strain lineages (`MYKROBE_LINEAGE.csv.gz`), train on everything else and test on that lineage. This checks the model isn't relying on strain markers.
   - **Benchmark:** every test sample is also called by the WHO mutation catalogue, 2nd edition (2023, v2.1, from [oxfordmmm/tuberculosis_amr_catalogues](https://github.com/oxfordmmm/tuberculosis_amr_catalogues), applied with [piezo](https://github.com/oxfordmmm/piezo)), and by the CRyPTIC catalogue. TB-Profiler bases its resistance calls on the WHO catalogue, so this is the catalogue-level equivalent of a TB-Profiler benchmark. A real TB-Profiler run would need the raw sequencing reads.

## Resistance level (MIC)

TB-Profiler and the WHO catalogue say resistant or susceptible. `mic.py` also asks how resistant, using the MICs measured on the CRyPTIC plates.

- **Isoniazid:** a ridge regression predicts the MIC and labels resistance as low-level (up to 2 doublings above the cut-off, typically *inhA*) or high-level (typically *katG*). On resistant samples from unseen countries it is within one doubling of the lab 87% of the time, against 79% for a yes/no call.
- **Other drugs:** the model did not beat a yes/no call, so the site shows the typical MIC for each resistance mutation instead, measured in samples where it is the only resistance mutation.

## Retraining in the cloud

Training runs on GitHub's machines through GitHub Actions. Nothing needs to be installed on your computer.

1. Open the repository on GitHub and go to the **Actions** tab.
2. Choose **Train and validate TB resistance model** and click **Run workflow**.
3. After about 10 minutes the results appear on the run's summary page and are saved to `results/`. The web page then updates with the new model.

The workflow also runs on its own whenever `train.py` or `requirements.txt` changes.

## Running locally (optional)

```bash
pip install -r requirements.txt
python train.py --data data --out results
python evaluate.py --data data --out results
python mic.py --data data --out results
python pick_examples.py --data data --out results
```

The data files must be in `data/`. Download links are in `.github/workflows/train.yml`.

## Known limitations

- Bedaquiline, linezolid, delamanid and clofazimine have fewer than 200 resistant samples, so those models are unreliable.
- Pyrazinamide is not on the UKMYC plates, so there is no pyrazinamide model.
- Mutations are only considered in genes the CRyPTIC catalogue links to each drug.
- The WHO catalogue was built partly from the same datasets used for testing here, so its benchmark numbers are optimistic.
- Mutations carried by only part of the bacterial population (mixed calls) are dropped, so heteroresistance can be missed.
