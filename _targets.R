library(targets)

tar_option_set(
  packages = c(
    "dplyr",
    "readr",
    "tibble",
    "tidyr",
    "purrr",
    "scales",
    "stringr",
    "ggplot2",
    "yaml",
    "jsonlite",
    "MOFA2"
  ),
  format = "rds"
)

invisible(lapply(list.files("R", pattern = "[.][Rr]$", full.names = TRUE), source))

list(
  tar_target(config, read_config("config/analysis.yaml"), cue = tar_cue(mode = "always")),
  tar_target(raw_multiomics, load_multiomics_data(config)),
  tar_target(qc_outputs, run_qc_workflow(raw_multiomics, config), format = "file"),
  tar_target(preprocessed_multiomics, preprocess_multiomics(raw_multiomics, config)),
  tar_target(eda_outputs, run_eda_workflow(preprocessed_multiomics, config), format = "file"),
  tar_target(mofa_run_full, train_mofa_run(preprocessed_multiomics, config, run_id = "run_01_full"), format = "file"),
  tar_target(mofa_run_filtered, train_mofa_run(preprocessed_multiomics, config, run_id = "run_02_filtered_proteins"), format = "file"),
  tar_target(downstream_outputs, run_downstream_workflow(preprocessed_multiomics, config, mofa_run_full, mofa_run_filtered), format = "file"),
  tar_target(evaluation_outputs, run_evaluation_workflow(preprocessed_multiomics, config), format = "file")
)
