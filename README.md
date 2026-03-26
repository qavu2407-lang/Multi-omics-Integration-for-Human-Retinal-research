# Multi-omics-Integration-for-Human-Retinal-research

A research-oriented hybrid workflow for retinal multi-omics studies spanning preprocessing, quality control, exploratory analysis, and downstream interpretation of externally trained MOFA+ models.

## Current architecture

The supported flow is intentionally split:

- `src/r/io.R`
- `src/r/preprocessing.R`
- `src/r/visualization.R`
- `src/r/mofa_01_io.R`
- `src/r/mofa_02_downstream.R`
- `_targets.R`
- `config/analysis.yaml`

Python notebooks are the source of truth for:

- MOFA+ model training
- robustness and evaluation analyses

R is the source of truth for:

- data loading
- QC summaries
- EDA and PCA plots
- downstream analysis of saved `.hdf5` MOFA models

## Repository layout

```text
.
├── README.md
├── config/
├── data/
├── notebooks/
│   └── python/
├── plots/
├── results/
├── scripts/
│   └── r/
└── src/
    └── r/
```

## Expected input format

Place the three CSV files in `data/`:

- `data/proteins.csv`
- `data/lipids.csv`
- `data/metabolites.csv`

Each file is expected to contain:

- one identifier column named `sample`
- one special row where `sample == "label"` containing biological group labels for each sample column
- all remaining rows representing molecular features

## Run the R pipeline

```r
Rscript scripts/r/run_pipeline.R
```

Phase-specific entrypoints:

```r
Rscript scripts/r/01_preprocess.R
Rscript scripts/r/02_train_and_downstream.R
Rscript scripts/r/03_evaluation.R
```

`scripts/r/03_evaluation.R` is an informational stub in the hybrid workflow because evaluation is handled in Python.

## External MOFA models

The R downstream step expects Python-trained MOFA models at the paths configured in `config/analysis.yaml`:

- `results/mofa/mofa/run_01_only_ARD.hdf5`
- `results/mofa/mofa/run_02_filtered_proteins.hdf5`

Or interactively:

```r
library(targets)
tar_make()
tar_visnetwork()
```

If you use the `targets` entrypoint, make sure the `targets` package is installed in your R environment.
