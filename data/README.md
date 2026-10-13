# data/

Raw and participant-level data files are not included in this public repository.

To reproduce the analyses, place approved local copies of the required input files in this directory. 

## Data availability

Data from the Kaiser Permanente Research Bank (KPRB), including the NMR metabolomics and linked electronic health record data used in this study, are available to qualified researchers through the [KPRB application process](https://researchbank.kaiserpermanente.org/for-researchers/apply-for-access/), subject to review and approval by KPRB and applicable institutional and data-use requirements.

## Expected files

All tables are keyed by `StudyID`. Files marked "created by" are written by an earlier notebook in this repository; the rest are KPRB extracts supplied by authorized users.

| File | Contents | Source |
|---|---|---|
| `assay.csv` | Needle-to-freezer-corrected NMR biomarker matrix, one row per participant | Created by `01_needle_to_freezer_correction` |
| `rowdata.csv` | Annotation for each NMR biomarker: name, unit, group, subgroup | Created by `01_needle_to_freezer_correction` |
| `meta.csv` | Participant metadata from the maplet colData (see [`preprocessing_inputs/README.md`](../preprocessing_inputs/README.md)) | Created by `01_needle_to_freezer_correction` |
| `meta_clin.csv` | `meta.csv` plus derived clinical predictors (labs, vitals, medications, smoking, PREVENT and CHARGE-AF ingredients, A1C variability, eGFR slope, fibrosis scores) used by all models | Created by `01_clin_preparation` |
| `a1c_history.csv` | Longitudinal HbA1c results | Created by `01_clin_preparation` |
| `egfr_history.csv` | Longitudinal creatinine-based eGFR | Created by `01_clin_preparation` |
| `fib4_history.csv` | Longitudinal FIB-4 with its AST, ALT, and platelet components | Created by `01_clin_preparation` |
| `labs.csv` | Laboratory results (test type, value, result date) | KPRB extract |
| `rx.csv` | Medication dispensings (generic name, dispensing date) | KPRB extract |
| `smoking.csv` | Smoking status | KPRB extract |
| `bmi.csv` | BMI measurements with dates | KPRB extract |
| `height.csv` | Height measurements with dates | KPRB extract |
| `weight.csv` | Weight measurements with dates | KPRB extract |
| `sbp.csv` | Systolic blood pressure measurements with dates | KPRB extract |
| `dbp.csv` | Diastolic blood pressure measurements with dates | KPRB extract |
| `prs.csv` | Polygenic scores per participant and PGS Catalog score | KPRB extract |
| `prs_catalog.csv` | Requested PGS Catalog IDs, each tagged with its target endpoint | Study PRS request |
| `healthcare_engagement/engage_list.rds` | Participant lists for three healthcare-engagement definitions (vitals only, vitals plus prescriptions, any contact) | Created by `01_healthcare_engagement`; used by `02_primary_incident_engagement` |
