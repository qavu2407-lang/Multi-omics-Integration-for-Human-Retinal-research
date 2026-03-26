normalise_mofa_columns <- function(tbl) {
  names(tbl) <- tolower(names(tbl))
  tbl
}

sample_id_key <- function(x) {
  gsub("^sample_", "", as.character(x))
}

extract_factors_table <- function(model, factors = NULL) {
  factor_tbl <- MOFA2::get_factors(model, factors = factors, as.data.frame = TRUE) %>%
    normalise_mofa_columns()

  if (all(c("sample", "factor", "value") %in% names(factor_tbl))) {
    id_cols <- intersect(c("sample", "group"), names(factor_tbl))
    factor_tbl <- factor_tbl %>%
      tidyr::pivot_wider(
        id_cols = dplyr::all_of(id_cols),
        names_from = factor,
        values_from = value
      ) %>%
      normalise_mofa_columns()
  }

  if ("sample" %in% names(factor_tbl)) {
    names(factor_tbl)[names(factor_tbl) == "sample"] <- "sample_id"
  }

  factor_tbl
}

summarise_variance_explained <- function(model, model_id) {
  variance_tbl <- MOFA2::get_variance_explained(model, as.data.frame = TRUE)$r2_per_factor %>%
    normalise_mofa_columns()

  variance_tbl %>%
    dplyr::group_by(view) %>%
    dplyr::summarise(r2 = sum(value), .groups = "drop") %>%
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
  factor_x_query <- factor_x
  factor_y_query <- factor_y
  factor_x <- tolower(factor_x_query)
  factor_y <- tolower(factor_y_query)

  plot_tbl <- extract_factors_table(model, factors = c(factor_x_query, factor_y_query)) %>%
    dplyr::mutate(sample_key = sample_id_key(sample_id)) %>%
    dplyr::left_join(
      sample_metadata %>% dplyr::mutate(sample_key = sample_id_key(sample_id)),
      by = "sample_key",
      suffix = c("", "_meta")
    ) %>%
    dplyr::mutate(
      label = dplyr::coalesce(group_meta, group),
      display_sample = paste0("S", sample_key)
    ) %>%
    dplyr::mutate(sample_id = dplyr::coalesce(sample_id, sample_id_meta)) %>%
    dplyr::select(-dplyr::any_of(c("sample_id_meta", "sample_key", "group", "group_meta")))

  ellipse_tbl <- plot_tbl %>%
    dplyr::filter(!is.na(label)) %>%
    dplyr::count(label, name = "n_samples") %>%
    dplyr::filter(n_samples >= 3)

  plot_obj <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes(x = .data[[factor_x]], y = .data[[factor_y]], color = label)
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = display_sample), size = 3, nudge_y = 0.05, show.legend = FALSE) +
    ggplot2::labs(
      title = paste0(run_label, ": ", factor_x_query, " vs ", factor_y_query),
      subtitle = "with 95% confidence ellipses",
      x = factor_x_query,
      y = factor_y_query,
      color = "Label"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  if (nrow(ellipse_tbl) > 0) {
    plot_obj <- plot_obj +
      ggplot2::stat_ellipse(
        data = dplyr::semi_join(plot_tbl, ellipse_tbl, by = "label"),
        ggplot2::aes(x = .data[[factor_x]], y = .data[[factor_y]], color = label),
        type = "norm",
        linewidth = 0.7,
        linetype = 2,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }

  plot_obj
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
  full_factor_plot <- plot_factor_scatter(
    full_model,
    preprocessed_multiomics$sample_metadata,
    "MOFA+ Run 1 (Full Data)"
  )
  filtered_factor_plot <- plot_factor_scatter(
    filtered_model,
    preprocessed_multiomics$sample_metadata,
    "MOFA+ Run 2 (Filtered Proteins)"
  )

  c(
    comparison_csv,
    save_plot(variance_plot, file.path(plot_dir, "variance_explained_comparison.png"), width = 8, height = 5),
    save_plot(full_weight_plot, file.path(plot_dir, "run_01_factor1_weights.png"), width = 9, height = 5),
    save_plot(filtered_weight_plot, file.path(plot_dir, "run_02_factor1_weights.png"), width = 9, height = 5),
    save_plot(full_factor_plot, file.path(plot_dir, "run_01_factors_F1_F2.png"), width = 8, height = 6),
    save_plot(filtered_factor_plot, file.path(plot_dir, "run_02_factors_F1_F2.png"), width = 8, height = 6)
  )
}
