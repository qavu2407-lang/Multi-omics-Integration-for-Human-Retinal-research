compute_missingness <- function(mat) {
  tibble::tibble(
    feature_id = rownames(mat),
    missing_fraction = rowMeans(is.na(mat)),
    variance = apply(mat, 1, stats::var, na.rm = TRUE)
  )
}

filter_matrix_by_missingness <- function(mat, feature_threshold) {
  feature_keep <- rowMeans(is.na(mat)) <= feature_threshold
  # Keep all samples for MOFA and only remove heavily missing features.
  mat[feature_keep, , drop = FALSE]
}

build_qc_summary <- function(multiomics, config) {
  purrr::imap_dfr(multiomics$views, function(mat, view_name) {
    tibble::tibble(
      view = view_name,
      n_features = nrow(mat),
      n_samples = ncol(mat),
      missing_values = sum(is.na(mat)),
      duplicated_features = anyDuplicated(rownames(mat)) > 0,
      duplicated_samples = anyDuplicated(colnames(mat)) > 0,
      median_feature_missingness = stats::median(rowMeans(is.na(mat))),
      median_sample_missingness = stats::median(colMeans(is.na(mat)))
    )
  })
}

run_qc_workflow <- function(multiomics, config) {
  qc_dir <- config$paths$qc_dir
  ensure_dir(qc_dir)

  summary_tbl <- build_qc_summary(multiomics, config)
  summary_path <- file.path(qc_dir, "qc_summary.csv")
  readr::write_csv(summary_tbl, summary_path)

  metadata_path <- file.path(qc_dir, "sample_metadata.csv")
  readr::write_csv(multiomics$sample_metadata, metadata_path)

  feature_missing_files <- purrr::imap_chr(multiomics$views, function(mat, view_name) {
    out_path <- file.path(qc_dir, paste0(view_name, "_feature_missingness.csv"))
    readr::write_csv(compute_missingness(mat), out_path)
    out_path
  })

  session_path <- write_session_metadata(file.path(qc_dir, "session_info.json"))
  c(summary_path, metadata_path, unname(feature_missing_files), session_path)
}

apply_variance_filter <- function(views, config) {
  vf <- config$data$variance_filter
  if (!isTRUE(vf$enabled)) {
    return(views)
  }

  target_view <- vf$view
  if (!target_view %in% names(views)) {
    stop(sprintf("Variance filter target view '%s' not found.", target_view))
  }

  mat <- views[[target_view]]
  ranked <- sort(apply(mat, 1, stats::var, na.rm = TRUE), decreasing = TRUE)
  keep_ids <- names(utils::head(ranked, min(vf$top_n, length(ranked))))
  views[[target_view]] <- mat[keep_ids, , drop = FALSE]
  views
}

preprocess_multiomics <- function(multiomics, config) {
  base_views <- purrr::map(multiomics$views, function(mat) {
    filter_matrix_by_missingness(
      mat,
      feature_threshold = config$data$feature_missingness_threshold
    )
  })

  filtered_views <- apply_variance_filter(base_views, config)

  sample_metadata <- multiomics$sample_metadata %>%
    dplyr::filter(sample_id %in% unique(unlist(purrr::map(base_views, colnames)))) %>%
    dplyr::arrange(sample_id)

  list(
    views = base_views,
    filtered_views = filtered_views,
    sample_metadata = sample_metadata,
    config = config
  )
}
