#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(targets))
tar_make(names = evaluation_outputs)
