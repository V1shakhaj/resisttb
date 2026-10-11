# TB drug-resistance model: results

Trained on the CRyPTIC dataset (lab-tested TB isolates with whole-genome sequences).
Each drug has its own model. 80% of samples were used for training and 20% were held back for testing.
Numbers below are on the held-back 20% only.

- **Sensitivity**: of the truly resistant samples, how many the method caught. Missing these is the dangerous error.
- **Specificity**: of the truly susceptible samples, how many were correctly called susceptible.
- **Catalogue**: the rule-based approach (known resistance mutations), with unknown mutations treated as susceptible.

| Drug | Samples (resistant) | Model sensitivity | Model specificity | Catalogue sensitivity | Catalogue specificity |
|---|---|---|---|---|---|
| Rifampicin | 10,396 (4,334) | 95.0% | 97.2% | 95.8% | 95.2% |
| Isoniazid | 10,432 (4,889) | 93.1% | 98.5% | 93.6% | 98.0% |
| Ethambutol | 7,755 (1,817) | 90.1% | 95.0% | 83.7% | 93.7% |
| Levofloxacin | 9,174 (1,615) | 82.0% | 97.2% | 81.4% | 97.3% |
| Moxifloxacin | 8,184 (1,336) | 88.8% | 94.2% | 83.1% | 94.5% |
| Amikacin | 10,414 (764) | 79.1% | 96.4% | 75.8% | 99.1% |
| Kanamycin | 10,793 (990) | 79.8% | 96.8% | 69.7% | 99.2% |
| Ethionamide | 9,508 (1,532) | 77.5% | 97.1% | 59.5% | 98.8% |
| Rifabutin | 11,500 (4,069) | 95.0% | 95.6% | no rules | no rules |
| Bedaquiline ⚠️ | 9,982 (94) | 68.4% | 81.0% | no rules | no rules |
| Linezolid ⚠️ | 8,575 (93) | 57.9% | 94.6% | no rules | no rules |
| Delamanid ⚠️ | 9,496 (99) | 50.0% | 64.0% | no rules | no rules |
| Clofazimine ⚠️ | 9,198 (139) | 64.3% | 79.9% | no rules | no rules |

⚠️ Fewer than 200 resistant samples to learn from. Treat these results as unreliable until more data is added.

"No rules" means the catalogue used in this dataset makes no resistance calls for that drug, so the model is the only predictor.

## Strongest resistance signals the model learned

- **Rifampicin**: `rpoB H445D`, `rpoB S450L`, `rpoB D435V`, `rpoB H445Y`, `rpoB S450W`
- **Isoniazid**: `katG S315T`, `katG S315N`, `ahpC c-54t`, `ahpC c-57t`, `fabG1 c-15t`
- **Ethambutol**: `embB M306V`, `embB E405D`, `embB Y319S`, `embB Q497R`, `embB Q445R`
- **Levofloxacin**: `gyrA D94Y`, `gyrA D94H`, `gyrB D461N`, `gyrA G88C`, `gyrA D94G`
- **Moxifloxacin**: `gyrA D94Y`, `gyrA D94G`, `gyrA D94N`, `gyrA G88C`, `gyrA D94H`
- **Amikacin**: `rrs a1401g`, `gid P75R`, `eis c-14t`, `gid A80P`, `gid L26F`
- **Kanamycin**: `rrs a1401g`, `eis c-14t`, `eis g-37t`, `eis g-10a`, `rrs c1402t`
- **Ethionamide**: `fabG1 c-15t`, `fabG1 L203L`, `ethA L35R`, `fabG1 c-34g`, `fabG1 V221V`
- **Rifabutin**: `rpoB H445D`, `rpoB S450L`, `rpoB H445Y`, `rpoB S450W`, `rpoB S450F`
- **Bedaquiline**: `atpE E61D`, `Rv0678 141_indel`, `Rv0678 M139I`, `mmpL3 L528L`, `Rv0678 Q51R`
- **Linezolid**: `rplC C154R`, `rrl g2814t`, `rrl g2270t`, `Rv1816 V57A`, `mmpL3 R923R`
- **Delamanid**: `ddn 362_indel`, `ddn L49P`, `ddn G34R`, `ddn Y65S`, `Rv1816 V57A`
- **Clofazimine**: `mmpL3 T284A`, `Rv3249c A17A`, `Rv0678 L117R`, `pepQ V45L`, `mmpL5 Y300!`

Research prototype. Not validated for clinical use.
