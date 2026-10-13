# NMR metabolomics for cardiovascular–kidney–metabolic risk stratification

Analysis code for evaluating whether NMR metabolomics improves risk prediction for cardiovascular–kidney–metabolic (CKM) disease in a real-world integrated healthcare system.

The study profiled 250 nuclear magnetic resonance (NMR) biomarkers in 54,933 participants of the Kaiser Permanente Research Bank (KPRB), linked to longitudinal electronic health records. For nine CKM endpoints, including type 2 diabetes, chronic kidney disease, MASLD, and cardiovascular outcomes, the notebooks:

- compare standard-of-care (SOC) Cox models with and without NMR biomarkers using out-of-fold cross-validation and paired C-index tests
- extend the comparison to disease-specific complications and cross-organ outcomes within T2D, CKD, and MASLD cohorts
- test robustness across demographic, anthropometric, and polygenic-risk strata and in complete-case, healthcare-engagement, and time-restriction sensitivity analyses
- assess calibration, high-risk enrichment, and decision-curve net benefit
- transport UK Biobank–trained models to KPRB
- illustrate NMR-enhanced risk estimates alongside individual patients' longitudinal clinical histories

This repository contains analysis code only. No participant data, derived results, or fitted models are included.

## Repository structure

| Path | Contents |
|---|---|
| `notebooks/` | Analysis notebooks (R kernel), numbered by stage |
| `helpers/` | R functions sourced by the notebooks |
| `data/` | Placeholder for analysis-ready KPRB inputs; see [`data/README.md`](data/README.md) |
| `preprocessing_inputs/` | Placeholder for raw NMR objects used by the preprocessing notebooks; see [`preprocessing_inputs/README.md`](preprocessing_inputs/README.md) |
| `outputs/` | Placeholder for generated results; see [`outputs/README.md`](outputs/README.md) |
| `models/` | Placeholder for UK Biobank–trained models used by the transport notebook; see [`models/README.md`](models/README.md) |
| `renv.lock`, `renv/`, `.Rprofile` | R package environment |

## Analysis workflow

| Stage | Notebook | Purpose |
|---|---|---|
| 01 Data preparation | `01_nmr_preprocessing` | NMR missingness filtering and imputation |
| | `01_needle_to_freezer_correction` | Correct NMR biomarkers for sample needle-to-freezer time; writes `data/assay.csv`, `data/meta.csv`, `data/rowdata.csv` |
| | `01_clin_preparation` | Build clinical predictors (labs, vitals, medications, smoking, PREVENT and CHARGE-AF ingredients, A1C variability, eGFR slope, fibrosis scores); writes `data/meta_clin.csv` and longitudinal histories |
| | `01_healthcare_engagement` | Define healthcare-engagement criteria (vitals, prescriptions, any contact) and write `data/healthcare_engagement/engage_list.rds` |
| 02 Primary endpoints | `02_primary_incident_prs_selection` | Select the endpoint-specific polygenic score used in the primary models |
| | `02_primary_incident_ckm` | Main out-of-fold Cox models (BASE, CLIN, SOC, each with and without PRS and NMR), C-index and delta C-index, and absolute risks |
| | `02_primary_incident_coef` | SOC+NMR model coefficients |
| | `02_primary_incident_prs_comparison` | PRS versus NMR C-index comparison |
| | `02_primary_calibration` | Calibration at 1, 3, and 5 years |
| | `02_primary_enrichment` | Observed incidence across predicted-risk percentiles |
| | `02_nmr_correlation` | NMR biomarker correlation matrix |
| | `02_primary_incident_missing_sens` | Complete-case sensitivity analysis |
| | `02_primary_incident_time_restrict` | Excluding events in the first 6 months or 1 year |
| | `02_primary_incident_engagement` | Healthcare-engagement sensitivity analysis |
| | `02_primary_ukbb_transport` | Apply UK Biobank–trained models to KPRB |
| 03 Subgroups | `03_primary_incident_subgroups` | C-index, enrichment, and calibration within demographic, anthropometric, and PRS strata |
| 04 Secondary outcomes | `04_secondary_diabetes`, `04_secondary_kidney`, `04_secondary_liver` | Complication and cross-organ endpoints within T2D, CKD, and MASLD cohorts |
| | `04_secondary_incident_coef`, `04_secondary_enrichment`, `04_secondary_calibration` | Coefficients, enrichment, and calibration for secondary endpoints |
| 05 Clinical utility | `05_dca_net_benefit` | Decision-curve net benefit |
| 06 Patient illustrations | `06_patient_plots` | Longitudinal clinical tracks with SOC and SOC+NMR risk estimates |

## Restore the R environment

The analyses were run with R 4.2.2. Install `renv`, then restore the package library from the lockfile:

```r
install.packages("renv")
renv::restore()
```

The project `.Rprofile` activates `renv` automatically when R is started in the repository root. The notebooks use the R Jupyter kernel from [IRkernel](https://github.com/IRkernel/IRkernel).

## Running the notebooks

1. Obtain approved data access (see [Data availability](#data-availability)) and place the input files listed in `data/README.md` and `preprocessing_inputs/README.md`.
2. Create the output folders:

   ```bash
   mkdir -p outputs/{figures,tables,main/absolute_risks} \
     outputs/secondary_{diabetes,kidney,liver}/absolute_risks \
     outputs/supp/{complete_case,engagement,time_restriction}
   ```

3. Run the stage 01 notebooks in order.
4. Run `02_primary_incident_prs_selection`, then `02_primary_incident_ckm`. The remaining stage 02, 03, 05, and 06 notebooks read the results of `02_primary_incident_ckm`, except the sensitivity notebooks and `02_nmr_correlation`, which only need stage 01.
5. Run the three stage 04 disease notebooks before the stage 04 coefficient, enrichment, and calibration notebooks.

The model-fitting notebooks (`02_primary_incident_ckm`, the three sensitivity notebooks, and the three stage 04 disease notebooks) fit many cross-validated Cox models and are the long-running steps.

`02_primary_incident_engagement` reads the engagement list written by `01_healthcare_engagement`. `02_primary_ukbb_transport` requires the UK Biobank–trained models, which are not distributed.

Notebooks are committed without execution outputs.

## Data availability

Data from the Kaiser Permanente Research Bank (KPRB), including the NMR metabolomics and linked electronic health record data used in this study, are available to qualified researchers through the [KPRB application process](https://researchbank.kaiserpermanente.org/for-researchers/apply-for-access/), subject to review and approval by KPRB and applicable institutional and data-use requirements.

UK Biobank data are available to approved researchers through the UK Biobank access process. This study used UK Biobank data under Application #588633.

No participant-level data, derived results, or fitted models are included in this repository.

## Contact

Jan Krumsiek — [jak2043@med.cornell.edu](mailto:jak2043@med.cornell.edu)
