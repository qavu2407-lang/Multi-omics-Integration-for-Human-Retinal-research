#!/usr/bin/env Rscript

load_pipeline_packages <- function() {
  pkgs <- c("magrittr", "dplyr", "readr", "tibble", "tidyr", "purrr", "scales", "stringr", "ggplot2", "yaml", "jsonlite", "MOFA2")
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

project_root <- function() {
  script_path <- NULL
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_arg <- args[grepl(file_arg, args)]
  if (length(script_arg) > 0) {
    script_path <- sub(file_arg, "", script_arg[[1]])
  }

  script_dir <- if (is.null(script_path)) "." else dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
  normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
}

source_r_modules <- function(root = project_root()) {
  for (path in list.files(file.path(root, "src", "r"), pattern = "[.][Rr]$", full.names = TRUE)) {
    source(path, local = .GlobalEnv)
  }
}

load_pipeline_context <- function(config_path = "config/analysis.yaml") {
  root <- project_root()
  setwd(root)
  load_pipeline_packages()
  source_r_modules(root)

  config <- read_config(config_path)
  raw_multiomics <- load_multiomics_data(config)
  preprocessed_multiomics <- preprocess_multiomics(raw_multiomics, config)

  list(
    root = root,
    config = config,
    raw_multiomics = raw_multiomics,
    preprocessed_multiomics = preprocessed_multiomics
  )
}
