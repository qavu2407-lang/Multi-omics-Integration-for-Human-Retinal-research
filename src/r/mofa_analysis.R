normalise_mofa_columns <- function(tbl) {
  names(tbl) <- tolower(names(tbl))
  tbl
}

prepare_mofa_object <- function(views, config) {
  model <- MOFA2::create_mofa(views)

  data_opts <- MOFA2::get_default_data_options(model)
  data_opts$scale_views <- isTRUE(config$mofa$scale_views)
  data_opts$scale_groups <- isTRUE(config$mofa$scale_groups)

  model_opts <- MOFA2::get_default_model_options(model)
  model_opts$num_factors <- config$mofa$factors
  model_opts$spikeslab_weights <- isTRUE(config$mofa$spikeslab_weights)
  model_opts$ard_factors <- isTRUE(config$mofa$ard_factors)
  model_opts$ard_weights <- isTRUE(config$mofa$ard_weights)

  train_opts <- MOFA2::get_default_training_options(model)
  train_opts$seed <- config$analysis$random_seed
  train_opts$maxiter <- config$mofa$maxiter
  train_opts$convergence_mode <- config$mofa$convergence_mode
  if ("dropR2" %in% names(train_opts)) {
    train_opts$dropR2 <- config$mofa$drop_r2
  }

  MOFA2::prepare_mofa(
    object = model,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
}

train_mofa_run <- function(preprocessed_multiomics, config, run_id = "run_01_full", seed = NULL, views_override = NULL) {
  views <- views_override %||% preprocessed_multiomics$views
  if (!is.null(seed)) {
    config$analysis$random_seed <- seed
  }

  mofa_dir <- config$paths$mofa_dir
  ensure_dir(mofa_dir)
  outfile <- file.path(mofa_dir, paste0(run_id, ".hdf5"))

  model <- prepare_mofa_object(views, config)
  MOFA2::run_mofa(model, outfile = outfile, use_basilisk = FALSE)
  outfile
}

load_mofa_model <- function(path) {
  MOFA2::load_model(path)
}

summarise_variance_explained <- function(model, model_id) {
  variance_tbl <- MOFA2::calculate_variance_explained(model)$r2_per_factor %>%
    normalise_mofa_columns()

  variance_tbl %>%
    dplyr::group_by(view) %>%
    dplyr::summarise(r2 = sum(r2), .groups = "drop") %>%
    dplyr::mutate(model = model_id)
}

plot_variance_explained_comparison <- function(model_full, model_filtered) {
  bind_tbl <- dplyr::bind_rows(
    summarise_variance_explained(model_full, "Full data"),
    summarise_variance_explained(model_filtered, "Filtered proteomics")
  )

  ggplot2::ggplot(bind_tbl, ggplot2::aes(x = view, y = r2, fill = model)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(
      title = "Variance explained by omics view",
      x = "Omics layer",
      y = "Total variance explained (R2)",
      fill = "MOFA run"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

extract_ranked_weights <- function(model, factor = "Factor1", top_n = 15) {
  weight_tbl <- MOFA2::get_weights(model, views = "all", factors = factor, as.data.frame = TRUE) %>%
    normalise_mofa_columns()
  value_col <- setdiff(names(weight_tbl), c("view", "feature", "factor"))[[1]]

  weight_tbl %>%
    dplyr::rename(weight = dplyr::all_of(value_col)) %>%
    dplyr::arrange(dplyr::desc(abs(weight))) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::slice_head(n = top_n)
}

plot_ranked_weights <- function(model, run_label, factor = "Factor1") {
  ranked <- extract_ranked_weights(model, factor = factor, top_n = 15)

  ggplot2::ggplot(ranked, ggplot2::aes(x = rank, y = weight, color = view, label = feature)) +
    ggplot2::geom_line(linewidth = 0.7, show.legend = FALSE) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_text(hjust = 0, nudge_x = 0.2, size = 3, show.legend = FALSE) +
    ggplot2::labs(
      title = paste(run_label, "-", factor, "weights"),
      x = "Absolute weight rank",
      y = "Weight"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::coord_cartesian(clip = "off")
}

plot_factor_scatter <- function(model, sample_metadata, run_label, factor_x = "Factor1", factor_y = "Factor2") {
  factors_tbl <- MOFA2::get_factors(model, factors = c(factor_x, factor_y), as.data.frame = TRUE) %>%
    normalise_mofa_columns()
  names(factors_tbl)[names(factors_tbl) == "sample"] <- "sample_id"

  factor_x <- tolower(factor_x)
  factor_y <- tolower(factor_y)

  plot_tbl <- factors_tbl %>%
    dplyr::left_join(sample_metadata, by = "sample_id")

  ggplot2::ggplot(plot_tbl, ggplot2::aes_string(x = factor_x, y = factor_y, color = "group", label = "sample_id")) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::stat_ellipse(type = "norm", linewidth = 0.7, linetype = 2, show.legend = FALSE) +
    ggplot2::geom_text(size = 3, nudge_y = 0.05, show.legend = FALSE) +
    ggplot2::labs(
      title = paste(run_label, "latent factor separation"),
      x = factor_x,
      y = factor_y,
      color = "Group"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

run_downstream_workflow <- function(preprocessed_multiomics, config, run_full_path, run_filtered_path) {
  plot_dir <- config$paths$plot_downstream_dir
  mofa_dir <- config$paths$mofa_dir
  ensure_dir(plot_dir)

  full_model <- load_mofa_model(run_full_path)
  filtered_model <- load_mofa_model(run_filtered_path)

  comparison_tbl <- dplyr::bind_rows(
    summarise_variance_explained(full_model, "run_01_full"),
    summarise_variance_explained(filtered_model, "run_02_filtered_proteins")
  )
  comparison_csv <- file.path(mofa_dir, "variance_explained_by_view.csv")
  readr::write_csv(comparison_tbl, comparison_csv)

  variance_plot <- plot_variance_explained_comparison(full_model, filtered_model)
  full_weight_plot <- plot_ranked_weights(full_model, "MOFA run 01")
  filtered_weight_plot <- plot_ranked_weights(filtered_model, "MOFA run 02")
  full_factor_plot <- plot_factor_scatter(full_model, preprocessed_multiomics$sample_metadata, "MOFA run 01")
  filtered_factor_plot <- plot_factor_scatter(filtered_model, preprocessed_multiomics$sample_metadata, "MOFA run 02")

  c(
    comparison_csv,
    save_plot(variance_plot, file.path(plot_dir, "variance_explained_comparison.png"), width = 8, height = 5),
    save_plot(full_weight_plot, file.path(plot_dir, "run_01_factor1_weights.png"), width = 9, height = 5),
    save_plot(filtered_weight_plot, file.path(plot_dir, "run_02_factor1_weights.png"), width = 9, height = 5),
    save_plot(full_factor_plot, file.path(plot_dir, "run_01_factors_F1_F2.png"), width = 8, height = 6),
    save_plot(filtered_factor_plot, file.path(plot_dir, "run_02_factors_F1_F2.png"), width = 8, height = 6)
  )
}
