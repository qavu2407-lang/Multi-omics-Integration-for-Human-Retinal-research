# Multi-omics-Integration-for-Human-Retinal-research

An end-to-end, research-oriented pipeline for multi-omics human retinal studies, spanning preprocessing, quality control, MOFA2 integration, downstream interpretation, and robustness assessment.

## R pipeline architecture

The original notebook workflow has been restructured into a modular R pipeline built with `targets` for reproducibility and cleaner separation of concerns:

- `R/io.R`: raw data ingestion, label extraction, and cross-view sample alignment
- `R/preprocessing.R`: missingness checks, filtering, imputation, and feature selection
- `R/visualization.R`: EDA, top-variable-feature plots, and PCA visualizations
- `R/mofa_analysis.R`: MOFA2 training and downstream latent-factor interpretation
- `R/evaluation.R`: repeated-seed stability and sample-size sensitivity analyses
- `_targets.R`: dependency-aware pipeline orchestration
- `config/analysis.yaml`: study configuration and analysis parameters

## Expected input format

Place the three CSV files in `data/`:

- `data/proteins.csv`
- `data/lipids.csv`
- `data/metabolites.csv`

Each file is expected to contain:

- one identifier column named `sample`
- one special row where `sample == "label"` containing biological group labels for each sample column
- all remaining rows representing molecular features

## Run the pipeline

```r
Rscript scripts/run_pipeline.R
```

Phase-specific entrypoints:

```r
Rscript scripts/01_preprocess.R
Rscript scripts/02_train_and_downstream.R
Rscript scripts/03_evaluation.R
```

Or interactively:

```r
library(targets)
tar_make()
tar_visnetwork()
```

## Research-facing improvements over the notebooks

- Centralized configuration instead of repeated hard-coded parameters
- Reusable functions instead of duplicated notebook cells
- Explicit quality-control outputs and session metadata
- Sample alignment checks across omics layers
- Defensible evaluation design using repeated-seed stability and sample-subset sensitivity
- A file-based pipeline better suited to collaboration, auditability, and manuscript preparation
