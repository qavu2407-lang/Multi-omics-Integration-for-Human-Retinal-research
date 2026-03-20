# Multi-omics-Integration-for-Human-Retinal-research

An end-to-end, research-oriented workflow for multi-omics studies spanning preprocessing, quality control, MOFA+ integration, downstream interpretation, and robustness assessment.

## Python pipeline architecture

The notebook logic now has a reusable Python package and CLI intended for running the same MOFA+ workflow on new datasets with a configuration file instead of ad hoc notebook edits:

- `src/python/mofapipeline/io.py`: data ingestion, label extraction, and cross-view sample alignment
- `src/python/mofapipeline/preprocessing.py`: QC summaries, missingness filtering, imputation, and variance filtering
- `src/python/mofapipeline/eda.py`: EDA plots, top-variable-feature visualizations, and PCA diagnostics
- `src/python/mofapipeline/mofa.py`: MOFA+ training and downstream interpretation
- `src/python/mofapipeline/evaluation.py`: repeated-seed stability and sample-size sensitivity analyses
- `src/python/mofapipeline/cli.py`: command-line entrypoint for staged or end-to-end execution
- `config/analysis_python.yaml`: reusable study configuration

## Legacy R pipeline

The original refactor into R remains available:

- `src/r/io.R`
- `src/r/preprocessing.R`
- `src/r/visualization.R`
- `src/r/mofa_analysis.R`
- `src/r/evaluation.R`
- `_targets.R`
- `config/analysis.yaml`

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
│   ├── python/
│   └── r/
├── src/
│   ├── python/
│   └── r/
└── tests/
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

## Install the Python pipeline

```bash
python3 -m pip install -e .[mofa]
```

## Run the Python pipeline

Full run:

```bash
python3 scripts/python/run_pipeline.py --config config/analysis_python.yaml --stage all
```

Or via the console entrypoint:

```bash
mofa-pipeline --config config/analysis_python.yaml --stage all
```

Run a single stage:

```bash
mofa-pipeline --config config/analysis_python.yaml --stage eda
mofa-pipeline --config config/analysis_python.yaml --stage train
mofa-pipeline --config config/analysis_python.yaml --stage evaluate
```

## Run the legacy R pipeline

```r
Rscript scripts/r/run_pipeline.R
```

Phase-specific entrypoints:

```r
Rscript scripts/r/01_preprocess.R
Rscript scripts/r/02_train_and_downstream.R
Rscript scripts/r/03_evaluation.R
```

Or interactively:

```r
library(targets)
tar_make()
tar_visnetwork()
```

## Research-facing improvements over the notebooks

- Centralized configuration instead of repeated hard-coded parameters
- Reusable stage-specific modules instead of duplicated notebook cells
- Explicit quality-control outputs and session metadata
- Sample alignment and label-consistency checks across omics layers
- Reproducible seeds and persistent MOFA artifacts for auditability
- Defensible evaluation design using repeated-seed stability and sample-subset sensitivity
- A package-based workflow better suited to collaboration, method reuse, and manuscript preparation
