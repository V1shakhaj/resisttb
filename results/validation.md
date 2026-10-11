# Validation: unseen countries, other studies, and the WHO catalogue

These tests check whether the AI model still works on data it was not built from,
and compare it with the rule-based catalogues on exactly the same samples.

- **Sensitivity (sens.)**: of the truly resistant samples, how many were caught. Missing these is the dangerous error.
- **Specificity (spec.)**: of the truly susceptible samples, how many were correctly called susceptible.
- **WHO catalogue**: WHO mutation catalogue, 2nd edition (2023, v2.1), applied with piezo. TB-Profiler bases its resistance calls on this catalogue, so this column is the catalogue-level stand-in for TB-Profiler.
- **CRyPTIC catalogue**: the catalogue shipped with the CRyPTIC data (v1.31).
- Catalogue mutations graded *uncertain* are counted as susceptible, as a clinical report would.

## Test 1: countries the model never saw

For each country, the model was trained on all other countries and then tested on that one.

| Drug | Samples (resistant) | AI model sens. | AI model spec. | WHO catalogue sens. | WHO catalogue spec. | CRyPTIC catalogue sens. | CRyPTIC catalogue spec. |
|---|---|---|---|---|---|---|---|
| Rifampicin | 10,396 (4,334) | 94.0% | 96.9% | 95.7% | 95.1% | 95.1% | 95.1% |
| Isoniazid | 10,432 (4,889) | 93.8% | 98.6% | 93.0% | 98.5% | 94.0% | 98.3% |
| Ethambutol | 7,755 (1,817) | 91.5% | 94.3% | 88.6% | 93.7% | 84.6% | 93.4% |
| Levofloxacin | 9,174 (1,615) | 83.5% | 97.1% | 84.0% | 97.1% | 81.9% | 97.3% |
| Moxifloxacin | 8,184 (1,336) | 86.3% | 93.9% | 86.1% | 93.9% | 80.5% | 94.4% |
| Amikacin | 10,414 (764) | 78.9% | 97.7% | 79.6% | 98.6% | 77.6% | 99.2% |
| Kanamycin | 10,793 (990) | 79.6% | 95.2% | 79.3% | 97.2% | 71.1% | 98.9% |
| Ethionamide | 9,508 (1,532) | 74.6% | 96.4% | 68.5% | 97.0% | 58.9% | 99.0% |
| Rifabutin | 11,500 (4,069) | 92.5% | 95.4% | no rules | no rules | no rules | no rules |
| Bedaquiline ⚠️ | 9,982 (94) | 51.1% | 83.1% | 5.3% | 99.8% | no rules | no rules |
| Linezolid ⚠️ | 8,575 (93) | 34.4% | 95.2% | 34.4% | 99.9% | no rules | no rules |
| Delamanid ⚠️ | 9,496 (99) | 56.6% | 61.2% | 11.1% | 100.0% | no rules | no rules |
| Clofazimine ⚠️ | 9,198 (139) | 41.0% | 87.6% | 2.2% | 99.8% | no rules | no rules |

### By country (at least 20 resistant and 20 susceptible samples)

| Country | Drug | Samples (resistant) | AI model sens. / spec. | WHO catalogue sens. / spec. |
|---|---|---|---|---|
| Brazil | Rifampicin | 275 (74) | 89.2% / 100.0% | 90.5% / 93.5% |
| China | Rifampicin | 990 (241) | 85.1% / 98.7% | 89.6% / 97.9% |
| Germany | Rifampicin | 837 (368) | 97.6% / 94.7% | 98.1% / 93.6% |
| India | Rifampicin | 1,319 (725) | 94.9% / 95.1% | 95.6% / 93.8% |
| Italy | Rifampicin | 1,523 (913) | 93.3% / 94.6% | 95.3% / 91.0% |
| Peru | Rifampicin | 2,297 (786) | 96.6% / 97.2% | 97.3% / 96.6% |
| South Africa | Rifampicin | 1,645 (627) | 91.9% / 98.3% | 95.2% / 95.1% |
| Sweden | Rifampicin | 354 (331) | 98.8% / 8.7% | 98.8% / 0.0% |
| Vietnam | Rifampicin | 1,071 (187) | 88.2% / 99.0% | 90.9% / 97.7% |
| Brazil | Isoniazid | 243 (82) | 87.8% / 98.8% | 86.6% / 99.4% |
| China | Isoniazid | 961 (270) | 84.1% / 98.6% | 79.6% / 99.0% |
| Germany | Isoniazid | 1,000 (470) | 97.2% / 99.4% | 97.2% / 99.4% |
| India | Isoniazid | 1,271 (720) | 94.7% / 96.6% | 94.6% / 96.6% |
| Italy | Isoniazid | 1,532 (933) | 91.9% / 98.5% | 91.9% / 98.2% |
| Peru | Isoniazid | 2,242 (933) | 96.0% / 99.2% | 95.7% / 98.9% |
| South Africa | Isoniazid | 1,699 (618) | 90.6% / 99.0% | 86.7% / 98.9% |
| Vietnam | Isoniazid | 1,035 (421) | 94.8% / 98.0% | 93.8% / 97.7% |
| Brazil | Ethambutol | 203 (23) | 95.7% / 99.4% | 95.7% / 99.4% |
| China | Ethambutol | 845 (146) | 82.9% / 99.1% | 83.6% / 99.3% |
| Germany | Ethambutol | 620 (84) | 96.4% / 89.2% | 94.0% / 88.4% |
| India | Ethambutol | 974 (461) | 93.7% / 91.2% | 84.4% / 91.2% |
| Italy | Ethambutol | 1,016 (216) | 94.0% / 86.9% | 86.1% / 87.4% |
| Peru | Ethambutol | 1,698 (328) | 88.1% / 96.2% | 92.1% / 95.8% |
| South Africa | Ethambutol | 1,456 (302) | 94.0% / 98.3% | 92.1% / 98.4% |
| Sweden | Ethambutol | 239 (170) | 90.0% / 63.8% | 92.4% / 14.5% |
| United States | Ethambutol | 57 (25) | 100.0% / 84.4% | 96.0% / 87.5% |
| Vietnam | Ethambutol | 647 (62) | 85.5% / 96.4% | 80.6% / 96.8% |
| China | Levofloxacin | 892 (117) | 67.5% / 99.5% | 67.5% / 99.5% |
| Germany | Levofloxacin | 863 (180) | 90.0% / 96.2% | 89.4% / 96.2% |
| India | Levofloxacin | 1,054 (385) | 89.4% / 92.8% | 90.1% / 92.4% |
| Italy | Levofloxacin | 1,397 (312) | 77.2% / 97.7% | 77.9% / 97.4% |
| Peru | Levofloxacin | 2,158 (137) | 65.7% / 99.0% | 65.7% / 99.2% |
| South Africa | Levofloxacin | 1,379 (231) | 87.0% / 99.1% | 88.3% / 98.9% |
| Sweden | Levofloxacin | 352 (219) | 95.0% / 42.1% | 95.9% / 42.1% |
| China | Moxifloxacin | 882 (105) | 70.5% / 99.0% | 71.4% / 98.6% |
| Germany | Moxifloxacin | 732 (121) | 89.3% / 91.2% | 87.6% / 91.3% |
| India | Moxifloxacin | 1,057 (355) | 94.6% / 87.5% | 93.0% / 87.3% |
| Italy | Moxifloxacin | 1,251 (181) | 80.1% / 90.0% | 80.1% / 89.8% |
| Peru | Moxifloxacin | 1,752 (127) | 69.3% / 98.3% | 69.3% / 98.5% |
| South Africa | Moxifloxacin | 1,247 (229) | 88.6% / 97.9% | 88.6% / 97.6% |
| Sweden | Moxifloxacin | 353 (189) | 94.7% / 35.4% | 96.8% / 36.0% |
| China | Amikacin | 986 (44) | 45.5% / 98.9% | 45.5% / 100.0% |
| Germany | Amikacin | 1,025 (45) | 75.6% / 94.0% | 75.6% / 97.1% |
| India | Amikacin | 1,339 (92) | 59.8% / 98.1% | 59.8% / 98.8% |
| Italy | Amikacin | 1,487 (116) | 80.2% / 97.2% | 79.3% / 98.9% |
| Peru | Amikacin | 2,399 (144) | 87.5% / 99.1% | 88.9% / 99.2% |
| South Africa | Amikacin | 1,494 (220) | 92.3% / 99.2% | 92.7% / 99.7% |
| Sweden | Amikacin | 354 (82) | 81.7% / 84.6% | 86.6% / 79.0% |

## Test 2: other studies

The model was trained on all CRyPTIC plate data and tested on samples from 7 other studies. These use different lab methods
(mostly MGIT liquid culture and solid media at WHO critical concentrations) and come from over 40 countries.
Patients who also appear in the CRyPTIC plate data were removed.

| Drug | Samples (resistant) | AI model sens. | AI model spec. | WHO catalogue sens. | WHO catalogue spec. | CRyPTIC catalogue sens. | CRyPTIC catalogue spec. |
|---|---|---|---|---|---|---|---|
| Rifampicin | 25,571 (6,403) | 89.8% | 99.3% | 93.0% | 98.7% | 92.4% | 98.7% |
| Isoniazid | 25,469 (7,916) | 91.7% | 98.9% | 90.2% | 99.1% | 92.4% | 98.8% |
| Ethambutol | 23,741 (3,161) | 81.6% | 95.5% | 82.6% | 95.5% | 84.2% | 95.3% |
| Levofloxacin | 536 (222) | 80.6% | 89.8% | 79.7% | 90.8% | 77.9% | 92.7% |
| Moxifloxacin | 6,053 (1,100) | 86.7% | 92.7% | 86.0% | 92.8% | 81.6% | 93.5% |
| Amikacin | 7,915 (659) | 69.0% | 99.1% | 74.7% | 98.6% | 69.5% | 99.6% |
| Kanamycin | 6,475 (594) | 66.7% | 97.8% | 71.7% | 98.0% | 59.3% | 99.3% |
| Ethionamide | 2,702 (785) | 53.4% | 87.7% | 48.2% | 88.9% | 38.5% | 90.6% |
| Rifabutin | 229 (131) | 96.9% | 83.7% | no rules | no rules | no rules | no rules |
| Bedaquiline ⚠️ | 96 (3) | 33.3% | 92.5% | 0.0% | 100.0% | no rules | no rules |
| Linezolid ⚠️ | 3,591 (46) | 76.1% | 87.7% | 69.6% | 99.5% | no rules | no rules |
| Delamanid ⚠️ | 97 (2) | 50.0% | 41.1% | 0.0% | 100.0% | no rules | no rules |
| Clofazimine ⚠️ | 2,245 (20) | 5.0% | 64.0% | 5.0% | 99.6% | no rules | no rules |

### By study (at least 20 resistant and 20 susceptible samples)

| Study | Drug | Samples (resistant) | AI model sens. / spec. | WHO catalogue sens. / spec. |
|---|---|---|---|---|
| CLIRES2020 | Rifampicin | 1,094 (510) | 72.7% / 95.4% | 73.5% / 94.7% |
| LID2015 | Rifampicin | 3,337 (395) | 90.1% / 99.4% | 93.9% / 98.5% |
| MUMBAI2018 | Rifampicin | 169 (123) | 98.4% / 87.0% | 99.2% / 84.8% |
| NEJM2018 | Rifampicin | 9,497 (2,395) | 91.2% / 99.4% | 95.1% / 98.9% |
| PHE2019 | Rifampicin | 3,977 (123) | 91.9% / 99.9% | 95.9% / 99.5% |
| SEQTREAT2020 | Rifampicin | 5,740 (2,480) | 92.0% / 99.4% | 94.6% / 98.7% |
| WHO2019 | Rifampicin | 1,757 (377) | 85.9% / 99.1% | 91.2% / 98.3% |
| CLIRES2020 | Isoniazid | 1,094 (585) | 89.9% / 92.1% | 89.6% / 92.5% |
| LID2015 | Isoniazid | 3,368 (649) | 86.6% / 99.0% | 85.7% / 98.9% |
| MUMBAI2018 | Isoniazid | 215 (171) | 97.7% / 75.0% | 98.2% / 75.0% |
| NEJM2018 | Isoniazid | 9,371 (2,774) | 91.9% / 99.2% | 89.5% / 99.5% |
| PHE2019 | Isoniazid | 3,964 (565) | 94.5% / 99.5% | 93.6% / 99.8% |
| SEQTREAT2020 | Isoniazid | 5,700 (2,676) | 92.5% / 99.1% | 90.8% / 99.2% |
| WHO2019 | Isoniazid | 1,757 (496) | 90.5% / 98.5% | 90.3% / 98.7% |
| LID2015 | Ethambutol | 3,360 (239) | 82.0% / 96.4% | 80.8% / 96.4% |
| MUMBAI2018 | Ethambutol | 93 (43) | 97.7% / 82.0% | 95.3% / 84.0% |
| NEJM2018 | Ethambutol | 9,180 (1,301) | 86.9% / 95.0% | 84.8% / 95.1% |
| PHE2019 | Ethambutol | 3,952 (101) | 77.2% / 99.4% | 84.2% / 99.5% |
| SEQTREAT2020 | Ethambutol | 5,399 (1,229) | 77.1% / 92.6% | 82.9% / 92.5% |
| WHO2019 | Ethambutol | 1,757 (248) | 74.2% / 95.0% | 68.5% / 94.6% |
| SEQTREAT2020 | Levofloxacin | 436 (160) | 89.4% / 93.1% | 90.0% / 93.8% |
| CLIRES2020 | Moxifloxacin | 2,297 (657) | 91.6% / 83.9% | 91.3% / 84.1% |
| LID2015 | Moxifloxacin | 582 (22) | 81.8% / 98.8% | 81.8% / 99.1% |
| NEJM2018 | Moxifloxacin | 868 (215) | 73.0% / 96.9% | 70.7% / 97.4% |
| PHE2019 | Moxifloxacin | 192 (37) | 81.1% / 96.1% | 81.1% / 94.2% |
| SEQTREAT2020 | Moxifloxacin | 684 (152) | 86.2% / 91.0% | 86.2% / 91.0% |
| CLIRES2020 | Amikacin | 2,166 (130) | 58.5% / 98.4% | 59.2% / 99.1% |
| LID2015 | Amikacin | 770 (61) | 63.9% / 99.9% | 65.6% / 98.9% |
| NEJM2018 | Amikacin | 2,031 (171) | 63.2% / 99.5% | 69.6% / 98.3% |
| SEQTREAT2020 | Amikacin | 1,923 (241) | 74.7% / 99.0% | 84.6% / 98.0% |
| WHO2019 | Amikacin | 840 (44) | 90.9% / 99.6% | 90.9% / 98.6% |

## Caveats

- **The WHO catalogue had a head start.** WHO built its 2023 catalogue from about 52,000 isolates, including the CRyPTIC data
  and several of the studies tested here. The catalogue has partly seen these samples before, so its numbers here are optimistic.
  The AI model never saw any test sample, so its numbers are the fairer estimate for new data.
- ⚠️ Fewer than 200 resistant samples in the training data. Treat as unreliable.
- The mutation table records insertions and deletions without their exact length, so WHO rules for frameshifts cannot fire.
  This mainly lowers WHO catalogue sensitivity for bedaquiline, clofazimine and delamanid.
- Only mutations in genes the CRyPTIC catalogue links to a drug are available. WHO rules for other genes are missed.
- A true TB-Profiler run needs raw sequencing reads, several terabytes for these samples, which this project does not process.
- Different studies use different lab methods and cut-offs, so some disagreement with the lab result is expected for any method.
- When the AI model and the WHO catalogue both disagree with the lab for the same group (for example Sweden),
  the cause is more likely that group's lab results (method, cut-off or sample selection) than either predictor.

Research prototype. Not validated for clinical use.
