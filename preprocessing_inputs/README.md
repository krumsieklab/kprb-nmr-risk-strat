# preprocessing_inputs/

Raw and participant-level data files are not included in this public repository.

To reproduce the analyses, place approved local copies of the required input files in this directory. Do not commit files from this directory unless they have been explicitly reviewed and approved for public release.

These files are only used by the two NMR preprocessing notebooks. Access follows the KPRB data availability statement in [`data/README.md`](../data/README.md).

## Expected files

| File | Contents | Used by |
|---|---|---|
| `kp55k_maplet.rds` | [maplet](https://github.com/krumsieklab/maplet) `SummarizedExperiment` holding the raw NMR biomarker matrix (assay), biomarker annotations (rowData), and participant metadata (colData), including demographics and the EHR-derived disease status and time variables for the primary and secondary endpoints | Input to `01_nmr_preprocessing` |
| `kp55k_maplet_preprocessed.rds` | The same maplet object after missingness filtering and imputation | Written by `01_nmr_preprocessing`; input to `01_needle_to_freezer_correction` |
| `kp55k_ntof.csv` | Needle-to-freezer time in hours for each sample | Input to `01_needle_to_freezer_correction` |
| `nmr_features.csv` | Annotation table for the NMR biomarkers: name, unit, group, subgroup | Reference only; not read by the notebooks |
