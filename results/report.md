# TB drug-resistance model: results

Trained on the CRyPTIC dataset (lab-tested TB isolates with whole-genome sequences).
Each drug has its own model. 80% of samples were used for training and 20% were held back for testing.
Numbers below are on the held-back 20% only.

- **Sensitivity**: of the truly resistant samples, how many the method caught. Missing these is the dangerous error.
- **Specificity**: of the truly susceptible samples, how many were correctly called susceptible.
- **Catalogue**: the rule-based approach (known resistance mutations), with unknown mutations treated as susceptible.

| Drug | Samples (resistant) | Model sensitivity | Model specificity | Catalogue sensitivity | Catalogue specificity |
|---|---|---|---|---|---|
| Rifampicin | 10,396 (4,334) | 95.0% | 97.0% | 95.8% | 95.2% |
| Isoniazid | 10,432 (4,889) | 92.5% | 98.4% | 92.8% | 98.2% |
| Ethambutol | 7,755 (1,817) | 90.6% | 94.1% | 77.4% | 94.4% |
| Levofloxacin | 9,174 (1,615) | 82.0% | 97.2% | 81.4% | 97.3% |
| Moxifloxacin | 8,184 (1,336) | 88.8% | 94.2% | 83.1% | 94.5% |
| Amikacin | 10,414 (764) | 77.1% | 97.8% | 75.8% | 99.1% |
| Kanamycin | 10,793 (990) | 79.8% | 98.1% | 69.7% | 99.2% |
| Ethionamide | 9,508 (1,532) | 76.8% | 97.1% | 59.5% | 98.8% |
| Rifabutin | 11,500 (4,069) | 94.6% | 95.8% | no rules | no rules |
| Bedaquiline | skipped: too little data (94 resistant samples, needs 200) | | | | |
| Linezolid | skipped: too little data (93 resistant samples, needs 200) | | | | |
| Delamanid | skipped: too little data (99 resistant samples, needs 200) | | | | |
| Clofazimine | skipped: too little data (139 resistant samples, needs 200) | | | | |

Drugs with fewer than 200 resistant samples get no model: with that little data the models learned noise.

"No rules" means the catalogue used in this dataset makes no resistance calls for that drug, so the model is the only predictor.

## Strongest resistance signals the model learned

- **Rifampicin**: `rpoB S450L`, `rpoB H445D`, `rpoB H445Y`, `rpoB D435V`, `rpoB S450W`
- **Isoniazid**: `katG S315T`, `katG S315N`, `fabG1 c-15t`, `katG D142N`, `fabG1 L203L`
- **Ethambutol**: `embB Q445R`, `embB M306V`, `embB Y319S`, `embB Q497R`, `embB Y334H`
- **Levofloxacin**: `gyrA D94Y`, `gyrA D94H`, `gyrB D461N`, `gyrA G88C`, `gyrA D94G`
- **Moxifloxacin**: `gyrA D94Y`, `gyrA D94G`, `gyrA D94N`, `gyrA G88C`, `gyrA D94H`
- **Amikacin**: `rrs a1401g`, `rrs g1484t`, `eis c-14t`, `rrs c1402t`, `rrs -26_indel`
- **Kanamycin**: `rrs a1401g`, `eis c-14t`, `eis g-37t`, `eis g-10a`, `rrs c1402t`
- **Ethionamide**: `fabG1 c-15t`, `fabG1 L203L`, `ethA L35R`, `fabG1 V221V`, `ethA 1055_indel`
- **Rifabutin**: `rpoB H445Y`, `rpoB S450L`, `rpoB H445D`, `rpoB S450W`, `rpoB Q432L`

Research prototype. Not validated for clinical use.
