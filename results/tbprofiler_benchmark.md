# End-to-end check with the real TB-Profiler

Raw sequencing reads for 81 patients (of 100 selected) were downloaded from ENA and run through
TB-Profiler. ResistTB's model was then applied to TB-Profiler's own mutation calls. None of these patients were
used to train ResistTB. Each cell is sensitivity / specificity against the lab result.

- **TB-Profiler**: resistant if TB-Profiler reports a resistance mutation for the drug.
- **ResistTB model**: resistant if ResistTB's model gives 50% or more.
- **TB-Profiler + ResistTB**: resistant if either says so. This is how the add-on works on the website.

| Drug | Samples (resistant) | TB-Profiler | ResistTB model | TB-Profiler + ResistTB |
|---|---|---|---|---|
| Rifampicin | 81 (33) | 97.0% / 97.9% | 93.9% / 97.9% | 97.0% / 97.9% |
| Isoniazid | 81 (45) | 97.8% / 100.0% | 84.4% / 100.0% | 97.8% / 100.0% |
| Ethambutol | 80 (28) | 100.0% / 86.5% | 96.4% / 86.5% | 100.0% / 86.5% |
| Levofloxacin | 1 (1) | 100.0% / n/a | 100.0% / n/a | 100.0% / n/a |
| Moxifloxacin | 31 (9) | 88.9% / 90.9% | 88.9% / 95.5% | 88.9% / 90.9% |
| Amikacin | 31 (4) | 100.0% / 100.0% | 100.0% / 70.4% | 100.0% / 70.4% |
| Kanamycin | 24 (6) | 100.0% / 100.0% | 100.0% / 100.0% | 100.0% / 100.0% |

## Caveats

- About 100 patients is a small sample: a difference of one or two patients moves these numbers by several percent.
- Patients were chosen to mix resistance profiles, so the share of resistant samples is higher than in most clinics.
- The lab results come from several studies and methods (mostly MGIT), which have their own error rates.
- TB-Profiler's database and ResistTB's model both draw partly on the WHO catalogue, which was built from some of these studies.

Research prototype. Not validated for clinical use.
