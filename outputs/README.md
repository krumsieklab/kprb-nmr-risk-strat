# outputs/

Generated analysis outputs are not included in this public repository.

Rerunning the notebooks may create spreadsheets, figures, logs, or other derived outputs in this directory. Review any generated file before deciding whether it is safe for public release.

## Layout created by the notebooks

Some notebooks write into these folders without creating them first. Create the tree once from the repository root before running:

```bash
mkdir -p outputs/{figures,tables,main/absolute_risks} \
  outputs/secondary_{diabetes,kidney,liver}/absolute_risks \
  outputs/supp/{complete_case,engagement,time_restriction}
```

| Folder | Contents |
|---|---|
| `main/` | Out-of-fold Cox results for the primary endpoints (`<endpoint>_<model>_oof.rds`) |
| `main/absolute_risks/` | Fold-specific baseline hazards and absolute risks for SOC and SOC+NMR models |
| `secondary_diabetes/`, `secondary_kidney/`, `secondary_liver/` | Out-of-fold results and absolute risks for complication endpoints within each disease cohort |
| `supp/complete_case/`, `supp/engagement/`, `supp/time_restriction/` | Sensitivity-analysis results |
| `c_index_summaries/` | C-index and delta C-index tables per endpoint |
| `tables/` | Summary tables (C-index statistics, coefficients, PRS selection, subgroup counts, calibration, enrichment) |
| `figures/` | Generated figures |

Model result files (`.rds`) contain participant-level out-of-fold predictions and must never be committed.
