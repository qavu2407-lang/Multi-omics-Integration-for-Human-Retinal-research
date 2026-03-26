save_plot <- function(plot_obj, path, width = 10, height = 6) {
  ensure_parent_dir(path)
  ggplot2::ggsave(filename = path, plot = plot_obj, width = width, height = height, dpi = 300)
  path
}

plot_distribution_panel <- function(mat, view_name) {
  values_tbl <- tibble::tibble(value = as.numeric(mat)) %>%
    dplyr::filter(!is.na(value))

  ggplot2::ggplot(values_tbl, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 50, fill = "#4C78A8", color = "white") +
    ggplot2::labs(
      title = paste("Global signal distribution:", stringr::str_to_title(view_name)),
      x = "Abundance / intensity",
      y = "Feature count"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

plot_top_variable_features <- function(mat, view_name, top_n = 20) {
  top_tbl <- tibble::tibble(
    feature_id = rownames(mat),
    variance = apply(mat, 1, stats::var, na.rm = TRUE)
  ) %>%
    dplyr::slice_max(order_by = variance, n = top_n) %>%
    dplyr::mutate(feature_id = factor(feature_id, levels = rev(feature_id)))

  ggplot2::ggplot(top_tbl, ggplot2::aes(x = variance, y = feature_id)) +
    ggplot2::geom_col(fill = "#F58518") +
    ggplot2::labs(
      title = paste("Top variable features:", stringr::str_to_title(view_name)),
      x = "Variance",
      y = "Feature"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

plot_layer_pca <- function(mat, sample_metadata, view_name, group_label = "group") {
  pca_input <- t(mat)
  pca_fit <- stats::prcomp(pca_input, center = TRUE, scale. = TRUE)

  pca_tbl <- tibble::as_tibble(pca_fit$x[, 1:2], rownames = "sample_id") %>%
    dplyr::left_join(sample_metadata, by = "sample_id")

  variance <- summary(pca_fit)$importance[2, 1:2] * 100

  ggplot2::ggplot(
    pca_tbl,
    ggplot2::aes(x = PC1, y = PC2, color = group)
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::stat_ellipse(type = "norm", linewidth = 0.7, linetype = 2, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = sample_id), size = 3, nudge_y = 0.05, show.legend = FALSE) +
    ggplot2::labs(
      title = paste("PCA:", stringr::str_to_title(view_name)),
      x = sprintf("PC1 (%.1f%%)", variance[[1]]),
      y = sprintf("PC2 (%.1f%%)", variance[[2]]),
      color = group_label
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

run_eda_workflow <- function(preprocessed_multiomics, config) {
  plot_dir <- config$paths$plot_eda_dir
  ensure_dir(plot_dir)

  outputs <- purrr::imap(preprocessed_multiomics$views, function(mat, view_name) {
    dist_plot <- plot_distribution_panel(mat, view_name)
    top_plot <- plot_top_variable_features(mat, view_name, config$analysis$top_variable_features)
    pca_plot <- plot_layer_pca(
      mat,
      preprocessed_multiomics$sample_metadata,
      view_name,
      group_label = config$data$group_variable %||% "group"
    )

    c(
      save_plot(dist_plot, file.path(plot_dir, paste0(view_name, "_distribution.png")), width = 8, height = 5),
      save_plot(top_plot, file.path(plot_dir, paste0(view_name, "_top_features.png")), width = 8, height = 5),
      save_plot(pca_plot, file.path(plot_dir, paste0("pca_", view_name, ".png")), width = 8, height = 6)
    )
  })

  unname(unlist(outputs))
}
