# MIC prediction: how resistant, not just whether

The MIC is the lowest drug concentration that stops the bacteria growing. Samples just above the cut-off
(low-level resistance) may still respond to a higher dose; samples far above it will not.

**Essential agreement**: the share of samples where the predicted MIC is within one doubling of the lab's.
This is the standard accuracy measure for MIC methods. Higher is better.

- **AI model**: ridge regression on mutations, one model per drug.
- **WHO yes/no call**: knowing only the WHO catalogue's resistant/susceptible call (what TB-Profiler reports),
  the best guess is the typical MIC for that call. The gap between this column and the AI model is the value
  of predicting the level, not just yes or no.
- **Most common MIC**: always guessing the drug's most common MIC. Any useful method must beat this.
- **Resistant samples** columns count only samples above the cut-off, where the level matters for treatment.

## Result

On resistant samples, the AI model clearly beats the yes/no call only for **isoniazid**, where it separates high-level
*katG* from low-level *inhA* resistance, and for kanamycin, which WHO no longer recommends. It roughly ties for most other
drugs and is worse for moxifloxacin. Isoniazid (✓) is therefore the only drug where the website predicts a level.
For other drugs, knowing which resistance mutation is present already says most of what there is to know about the
level, so the website shows the typical MIC range for each mutation instead. It is measured in samples where that is
the only resistance-linked mutation for the drug (at least 10 such samples).
Low-level resistance means up to 2 doublings above the cut-off; high-level means further above.

## Unseen countries (each country predicted by a model trained on the others)

| Drug | Samples | Cut-off (ECOFF, mg/L) | Resistant samples: AI model | Resistant samples: WHO yes/no | All samples: AI model | All samples: WHO yes/no | All samples: most common MIC |
|---|---|---|---|---|---|---|---|
| Rifampicin | 10,396 | 0.5 | 91.1% | 93.1% | 93.8% | 90.5% | 32.7% |
| **Isoniazid** ✓ | 10,977 | 0.1 | 86.7% | 79.1% | 92.7% | 88.9% | 36.0% |
| Ethambutol | 8,940 | 2 | 84.7% | 84.0% | 91.7% | 89.4% | 65.1% |
| Levofloxacin | 9,174 | 1 | 70.0% | 67.6% | 88.8% | 88.1% | 78.4% |
| Moxifloxacin | 8,184 | 1 | 61.6% | 86.1% | 78.3% | 81.1% | 69.7% |
| Amikacin | 10,414 | 1 | 77.1% | 77.4% | 95.9% | 95.2% | 90.3% |
| Kanamycin | 10,793 | 4 | 74.7% | 69.1% | 96.4% | 94.6% | 90.8% |
| Ethionamide | 10,261 | 2 | 64.1% | 61.8% | 88.3% | 87.0% | 75.4% |
| Rifabutin | 11,500 | 0.12 | 83.4% | no rules | 91.1% | no rules | 64.6% |
| Bedaquiline | 9,982 | 0.25 | 0.0% | 1.1% | 92.1% | 92.4% | 92.3% |
| Linezolid | 8,575 | 1 | 30.1% | 34.4% | 85.9% | 88.3% | 88.0% |
| Delamanid | 9,496 | 0.12 | 0.0% | 11.1% | 95.1% | 95.2% | 95.1% |
| Clofazimine | 9,198 | 0.25 | 0.0% | 0.0% | 93.6% | 93.7% | 93.7% |

## Random 80/20 split

| Drug | Samples | Cut-off (ECOFF, mg/L) | Resistant samples: AI model | Resistant samples: WHO yes/no | All samples: AI model | All samples: WHO yes/no | All samples: most common MIC |
|---|---|---|---|---|---|---|---|
| Rifampicin | 10,396 | 0.5 | 92.2% | 93.5% | 94.0% | 90.1% | 40.3% |
| **Isoniazid** ✓ | 10,977 | 0.1 | 89.0% | 78.4% | 94.2% | 88.5% | 41.5% |
| Ethambutol | 8,940 | 2 | 86.3% | 84.5% | 92.9% | 90.0% | 65.4% |
| Levofloxacin | 9,174 | 1 | 67.8% | 66.6% | 88.7% | 88.0% | 78.5% |
| Moxifloxacin | 8,184 | 1 | 73.4% | 89.1% | 82.3% | 81.0% | 68.8% |
| Amikacin | 10,414 | 1 | 77.1% | 74.5% | 95.6% | 94.4% | 90.1% |
| Kanamycin | 10,793 | 4 | 73.7% | 69.7% | 97.0% | 95.6% | 90.8% |
| Ethionamide | 10,261 | 2 | 67.4% | 60.6% | 88.8% | 86.5% | 75.0% |
| Rifabutin | 11,500 | 0.12 | 87.7% | no rules | 94.0% | no rules | 64.6% |
| Bedaquiline | 9,982 | 0.25 | 0.0% | 0.0% | 92.3% | 92.2% | 92.0% |
| Linezolid | 8,575 | 1 | 31.6% | 31.6% | 88.7% | 88.8% | 88.5% |
| Delamanid | 9,496 | 0.12 | 0.0% | 30.0% | 95.0% | 95.2% | 94.8% |
| Clofazimine | 9,198 | 0.25 | 3.6% | 0.0% | 94.5% | 94.2% | 94.2% |

## Caveats

- MICs were measured on CRyPTIC UKMYC plates. Other lab methods give somewhat different values.
- The two plate versions test different ranges, so values beyond the shared range are grouped in the end bins (≤ and >).
- Effects of mutations are treated as additive on the doubling scale, which is a simplification.
- Dosing decisions need a specialist and local guidance. Research prototype, not for clinical use.
