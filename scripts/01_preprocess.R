#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(targets))
tar_make(names = c(qc_outputs, eda_outputs))
