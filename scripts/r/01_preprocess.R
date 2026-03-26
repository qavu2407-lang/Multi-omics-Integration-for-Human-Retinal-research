#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else "scripts/r/01_preprocess.R"
source(file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)), "common.R"))

ctx <- load_pipeline_context()

message("Running preprocessing workflow...")
qc_outputs <- run_qc_workflow(ctx$raw_multiomics, ctx$config)
eda_outputs <- run_eda_workflow(ctx$preprocessed_multiomics, ctx$config)

message("QC outputs:")
message(paste(qc_outputs, collapse = "\n"))
message("EDA outputs:")
message(paste(eda_outputs, collapse = "\n"))
