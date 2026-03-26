#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else "scripts/r/02_train_and_downstream.R"
source(file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)), "common.R"))

ctx <- load_pipeline_context()

message("Running downstream analysis from externally trained MOFA models...")
model_paths <- expected_mofa_model_paths(ctx$config)
mofa_run_full <- require_existing_mofa_model(model_paths$run_full)
mofa_run_filtered <- require_existing_mofa_model(model_paths$run_filtered)
downstream_outputs <- run_downstream_workflow(
  ctx$preprocessed_multiomics,
  ctx$config,
  mofa_run_full,
  mofa_run_filtered
)

message("Loaded MOFA models:")
message(paste(c(mofa_run_full, mofa_run_filtered), collapse = "\n"))
message("Downstream outputs:")
message(paste(downstream_outputs, collapse = "\n"))
