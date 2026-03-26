#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else "scripts/r/03_evaluation.R"
source(file.path(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)), "common.R"))

ctx <- load_pipeline_context()

message("Evaluation is handled in the Python workflow for this hybrid pipeline.")
message("No R-side MOFA evaluation was run.")
