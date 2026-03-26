expected_mofa_model_paths <- function(config) {
  list(
    run_full = config$paths$run_full_model,
    run_filtered = config$paths$run_filtered_model
  )
}

require_existing_mofa_model <- function(path) {
  if (!file.exists(path)) {
    stop(
      sprintf(
        paste(
          "Expected an externally trained MOFA model at '%s'.",
          "Train the model in Python first, then rerun the R downstream script."
        ),
        path
      )
    )
  }

  path
}

load_mofa_model <- function(path) {
  suppressWarnings(MOFA2::load_model(path))
}
