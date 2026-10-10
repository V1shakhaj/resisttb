# ResistTB Predictor

Predicts drug-resistant tuberculosis from the mutations found when a TB sample is sequenced.
Research prototype. Not for clinical use.

- **Live demo:** https://v1shakhaj.github.io/resisttb/ (`index.html`, served by GitHub Pages)
- **Training code:** `train.py`
- **Latest results:** [`results/report.md`](results/report.md)

## How it works

1. **Data:** the [CRyPTIC](https://www.crypticproject.org/) June 2022 release: about 12,000 TB isolates with whole-genome sequences and lab-measured resistance to 13 drugs (UKMYC microtitre plates).
   - `EFFECTS.csv.gz`: each sample's mutations in the genes linked to each drug, plus the catalogue's call for each mutation.
   - `UKMYC_PHENOTYPES.csv.gz`: the lab result (resistant or susceptible) per sample and drug. Only high- and medium-quality readings are used.
2. **Model:** one logistic regression per drug with an L1 penalty, which keeps only the mutations that matter. Each mutation gets a weight, so every prediction can be explained.
3. **Evaluation:** a stratified 80/20 split. The 20% test set is never used for training. Results are compared with the catalogue (known resistance mutations, with unknown mutations counted as susceptible).

## Retraining in the cloud

Training runs on GitHub's machines through GitHub Actions. Nothing needs to be installed on your computer.

1. Open the repository on GitHub and go to the **Actions** tab.
2. Choose **Train TB resistance model** and click **Run workflow**.
3. After about 5 minutes the results appear on the run's summary page and are saved to `results/`. The web page then updates with the new model.

The workflow also runs on its own whenever `train.py` or `requirements.txt` changes.

## Running locally (optional)

```bash
pip install -r requirements.txt
python train.py --data data --out results
```

The two data files must be in `data/`. Download links are in `.github/workflows/train.yml`.

## Known limitations

- Bedaquiline, linezolid, delamanid and clofazimine have fewer than 200 resistant samples, so those models are unreliable.
- Pyrazinamide is not on the UKMYC plates, so there is no pyrazinamide model.
- Mutations are only considered in genes the CRyPTIC catalogue links to each drug.
- A random split can overestimate accuracy when closely related strains appear in both training and test sets. A next step is to split by country or by strain lineage.
