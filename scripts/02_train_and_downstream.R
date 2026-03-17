#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(targets))
tar_make(names = c(mofa_run_full, mofa_run_filtered, downstream_outputs))
