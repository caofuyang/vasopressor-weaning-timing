# Timing of the first vasopressor weaning attempt and mortality in sepsis

Analysis code for a target trial emulation in MIMIC-IV with external validation in eICU-CRD.

**Question.** Among adults with sepsis receiving norepinephrine, is attempting weaning within six hours of first meeting objective weaning criteria associated with 28-day mortality, compared with delaying the attempt?

**Design.** Target trial emulation. Eligibility (t0) = first hour with norepinephrine ≤0.10 µg/kg/min, MAP ≥65 mmHg in that hour and the next, and no dose reduction ≥10% during the preceding six hours. Landmark = t0 + 6 h. Exposure = norepinephrine reduction ≥25% or complete discontinuation between t0 and t0 + 6 h. Primary analysis = propensity-score IPTW with weighted Cox regression, multiple imputation (20 datasets, Rubin's rules), robust variance.

## Data access

The patient-level data are **not** included in this repository and cannot be redistributed.

- MIMIC-IV v3.1 — https://physionet.org/content/mimiciv/ (credentialed access, data use agreement required)
- eICU-CRD v2.0 — https://physionet.org/content/eicu-crd/ (credentialed access)

Both require completion of the CITI "Data or Specimens Only Research" training and signing of a data use agreement.

## Repository layout

```
sql/
  03_cohort_v2.sql      MIMIC-IV cohort: sepsis diagnosis codes, eligibility, exposure, landmark, covariates, export
  07_eicu_cohort.sql    eICU-CRD cohort built with identical definitions
  15_sepsis_filter.sql  Impact assessment of the ICD-based sepsis restriction (run first)
code/
  04_analysis.py        Primary MIMIC analysis (IPTW, weighted Cox, multiple imputation, sensitivity analyses)
  05_table1.py          Table 1 (baseline characteristics by exposure group)
  08_eicu_analysis.py   eICU-CRD analysis
  09_mimic_hosp.py      MIMIC-IV re-analysis with in-hospital mortality (comparability with eICU)
  10_figure2.py         Figure 2, IPTW-weighted Kaplan-Meier curves
  11_figure1.py         Figure 1, cohort flow diagram
  12_subgroups.py       Subgroup analyses with interaction tests
  13_tables.py          Tables 1-3 rendered in English
  20_key_numbers.py     Single source of truth: computes every reported number into results/key_numbers.json
  21_verify.py          Build-time consistency gate (numbers and section order)
```

## Reproduction

```bash
# 1. Build the cohorts (DuckDB 1.5.x)
duckdb tte.duckdb -c ".read sql/15_sepsis_filter.sql"   # optional: check the impact of the sepsis restriction
duckdb tte.duckdb -c ".read sql/03_cohort_v2.sql"       # writes analytic_v2.csv
duckdb eicu.duckdb -c ".read sql/07_eicu_cohort.sql"    # writes eicu_analytic.csv

# 2. Primary and secondary analyses
python code/20_key_numbers.py     # results/key_numbers.json
python code/04_analysis.py        # primary estimate, sensitivity analyses
python code/08_eicu_analysis.py   # external validation
python code/09_mimic_hosp.py      # in-hospital mortality, MIMIC

# 3. Tables, figures, consistency check
python code/05_table1.py
python code/12_subgroups.py
python code/13_tables.py
python code/10_figure2.py
python code/11_figure1.py
python code/21_verify.py          # must print "all consistent" before submission
```

Python 3.10; see `requirements.txt`.

## Key results

| Analysis | Cohort | Estimate |
|---|---|---|
| Primary (28-day mortality) | MIMIC-IV, n = 2,781 | HR 0.72 (95% CI 0.63–0.81) |
| In-hospital mortality | MIMIC-IV, n = 2,780 | HR 0.74 (0.64–0.84) |
| External validation | eICU-CRD, n = 532 | HR 0.85 (0.63–1.13) |
| E-value | — | 2.14 |

## Citation

Manuscript in preparation. Please cite this repository together with the paper once published.

## Licence

MIT (see `LICENSE`). The licence covers the code only, not the underlying data.
