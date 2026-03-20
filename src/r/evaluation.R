subset_views_by_samples <- function(views, sample_ids) {
  purrr::map(views, ~ .x[, sample_ids, drop = FALSE])
}

seed_stability_analysis <- function(preprocessed_multiomics, config) {
  n_runs <- config$evaluation$repeated_seed_runs
  threshold <- config$evaluation$stability_correlation_threshold
  base_seed <- config$analysis$random_seed

  run_paths <- purrr::map_chr(seq_len(n_runs), function(i) {
    train_mofa_run(
      preprocessed_multiomics = preprocessed_multiomics,
      config = config,
      run_id = sprintf("stability_seed_%02d", i),
      seed = base_seed + i - 1
    )
  })

  models <- purrr::map(run_paths, load_mofa_model)
  variance_tbl <- purrr::imap_dfr(models, function(model, idx) {
    summarise_variance_explained(model, sprintf("seed_%02d", idx))
  })

  factor_tbls <- purrr::map(models, ~ MOFA2::get_factors(.x, as.data.frame = TRUE) %>% normalise_mofa_columns())
  ref_factor_tbl <- factor_tbls[[1]]
  factor_cols <- grep("^factor", names(ref_factor_tbl), value = TRUE)

  replication_tbl <- purrr::map_dfr(factor_cols, function(factor_name) {
    ref_vec <- ref_factor_tbl[[factor_name]]
    corrs <- purrr::map_dbl(factor_tbls[-1], function(tbl) {
      corr_vals <- stats::cor(ref_vec, as.matrix(tbl[, factor_cols, drop = FALSE]), use = "pairwise.complete.obs")
      max(abs(corr_vals), na.rm = TRUE)
    })
    tibble::tibble(
      factor = factor_name,
      mean_abs_correlation = mean(corrs),
      replication_rate = mean(corrs >= threshold)
    )
  })

  list(
    run_paths = run_paths,
    variance = variance_tbl,
    replication = replication_tbl
  )
}

sample_size_sensitivity_analysis <- function(preprocessed_multiomics, config) {
  fractions <- config$evaluation$sample_size_grid
  n_reps <- config$evaluation$replicates_per_size
  sample_ids <- preprocessed_multiomics$sample_metadata$sample_id
  base_seed <- config$analysis$random_seed

  purrr::map_dfr(fractions, function(frac) {
    target_n <- max(3, floor(length(sample_ids) * frac))
    purrr::map_dfr(seq_len(n_reps), function(rep_id) {
      set.seed(base_seed + rep_id + target_n)
      chosen_ids <- sort(sample(sample_ids, size = target_n, replace = FALSE))
      subset_views <- subset_views_by_samples(preprocessed_multiomics$views, chosen_ids)
      subset_meta <- preprocessed_multiomics$sample_metadata %>%
        dplyr::filter(sample_id %in% chosen_ids) %>%
        dplyr::mutate(sample_id = factor(sample_id, levels = chosen_ids)) %>%
        dplyr::arrange(sample_id) %>%
        dplyr::mutate(sample_id = as.character(sample_id))

      subset_data <- list(
        views = subset_views,
        sample_metadata = subset_meta,
        config = config
      )

      run_path <- train_mofa_run(
        preprocessed_multiomics = subset_data,
        config = config,
        run_id = sprintf("subset_%02d_frac_%s_rep_%02d", target_n, gsub("[.]", "_", as.character(frac)), rep_id),
        seed = base_seed + rep_id
      )

      summarise_variance_explained(load_mofa_model(run_path), sprintf("n_%02d_rep_%02d", target_n, rep_id)) %>%
        dplyr::mutate(sample_fraction = frac, n_samples = target_n, replicate = rep_id)
    })
  })
}

plot_seed_stability <- function(seed_results) {
  ggplot2::ggplot(seed_results$variance, ggplot2::aes(x = model, y = r2, color = view, group = view)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      title = "Repeated-seed stability of variance explained",
      x = "Repeated MOFA fit",
      y = "Total variance explained (R2)",
      color = "View"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}

plot_factor_replication <- function(seed_results) {
  ggplot2::ggplot(seed_results$replication, ggplot2::aes(x = factor, y = replication_rate)) +
    ggplot2::geom_col(fill = "#54A24B") +
    ggplot2::geom_hline(yintercept = 0.7, linetype = 2, color = "#E45756") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      title = "Factor replication across repeated seeds",
      x = "Latent factor",
      y = "Replication rate"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

plot_sample_sensitivity <- function(sensitivity_tbl) {
  summary_tbl <- sensitivity_tbl %>%
    dplyr::group_by(view, sample_fraction, n_samples) %>%
    dplyr::summarise(mean_r2 = mean(r2), sd_r2 = stats::sd(r2), .groups = "drop")

  ggplot2::ggplot(summary_tbl, ggplot2::aes(x = n_samples, y = mean_r2, color = view)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = n_samples, ymin = mean_r2 - sd_r2, ymax = mean_r2 + sd_r2, fill = view),
      alpha = 0.15,
      linewidth = 0,
      color = NA,
      inherit.aes = FALSE,
      data = summary_tbl
    ) +
    ggplot2::labs(
      title = "Sample-size sensitivity analysis",
      x = "Number of samples",
      y = "Mean variance explained (R2)",
      color = "View",
      fill = "View"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

run_evaluation_workflow <- function(preprocessed_multiomics, config) {
  eval_dir <- config$paths$evaluation_dir
  plot_dir <- config$paths$plot_evaluation_dir
  ensure_dir(eval_dir)
  ensure_dir(plot_dir)

  seed_results <- seed_stability_analysis(preprocessed_multiomics, config)
  sensitivity_tbl <- sample_size_sensitivity_analysis(preprocessed_multiomics, config)

  seed_variance_csv <- file.path(eval_dir, "repeated_seed_variance.csv")
  seed_replication_csv <- file.path(eval_dir, "repeated_seed_replication.csv")
  sensitivity_csv <- file.path(eval_dir, "sample_size_sensitivity.csv")

  readr::write_csv(seed_results$variance, seed_variance_csv)
  readr::write_csv(seed_results$replication, seed_replication_csv)
  readr::write_csv(sensitivity_tbl, sensitivity_csv)

  stability_plot <- plot_seed_stability(seed_results)
  replication_plot <- plot_factor_replication(seed_results)
  sensitivity_plot <- plot_sample_sensitivity(sensitivity_tbl)

  c(
    seed_variance_csv,
    seed_replication_csv,
    sensitivity_csv,
    save_plot(stability_plot, file.path(plot_dir, "repeated_seed_variance.png"), width = 10, height = 6),
    save_plot(replication_plot, file.path(plot_dir, "factor_replication.png"), width = 8, height = 5),
    save_plot(sensitivity_plot, file.path(plot_dir, "sample_size_sensitivity.png"), width = 8, height = 5)
  )
}
