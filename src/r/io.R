read_omics_layer <- function(path, view_name, id_column = "sample", label_token = "label") {
  raw_tbl <- readr::read_csv(path, show_col_types = FALSE)

  if (!id_column %in% names(raw_tbl)) {
    stop(sprintf("Column '%s' not found in %s.", id_column, path))
  }

  label_row <- raw_tbl %>%
    dplyr::filter(.data[[id_column]] == label_token)

  if (nrow(label_row) != 1) {
    stop(sprintf("Expected exactly one label row in %s, found %s.", path, nrow(label_row)))
  }

  feature_tbl <- raw_tbl %>%
    dplyr::filter(.data[[id_column]] != label_token) %>%
    dplyr::mutate(dplyr::across(-dplyr::all_of(id_column), readr::parse_number))

  raw_sample_ids <- setdiff(names(feature_tbl), id_column)
  sample_ids <- paste0("sample_", raw_sample_ids)
  feature_ids <- feature_tbl[[id_column]]

  if (anyDuplicated(feature_ids)) {
    stop(sprintf("Duplicated feature identifiers detected in %s.", path))
  }

  mat <- feature_tbl %>%
    tibble::column_to_rownames(id_column) %>%
    as.matrix()
  storage.mode(mat) <- "double"
  colnames(mat) <- sample_ids

  sample_metadata <- tibble::tibble(
    sample_id = sample_ids,
    group = as.character(unlist(label_row[1, raw_sample_ids, drop = TRUE]))
  )

  list(
    view = view_name,
    matrix = mat,
    sample_metadata = sample_metadata,
    raw_dimensions = dim(raw_tbl)
  )
}

merge_sample_metadata <- function(layer_list) {
  combined_metadata <- purrr::map_dfr(layer_list, function(layer) {
    layer$sample_metadata %>%
      dplyr::mutate(source_view = layer$view)
  })

  conflicts <- combined_metadata %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::count(sample_id, name = "n_groups") %>%
    dplyr::filter(n_groups > 1)

  if (nrow(conflicts) > 0) {
    stop("Conflicting sample labels were found across omics layers.")
  }

  combined_metadata %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      group = dplyr::first(group),
      available_views = paste(sort(unique(source_view)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::arrange(sample_id)
}

collect_multiomics_layers <- function(layer_list) {
  list(
    views = purrr::set_names(purrr::map(layer_list, "matrix"), purrr::map_chr(layer_list, "view")),
    sample_metadata = merge_sample_metadata(layer_list)
  )
}

load_multiomics_data <- function(config) {
  paths <- config$paths
  settings <- config$data

  layers <- list(
    read_omics_layer(paths$proteomics, "proteomics", settings$id_column, settings$label_token),
    read_omics_layer(paths$lipidomics, "lipidomics", settings$id_column, settings$label_token),
    read_omics_layer(paths$metabolomics, "metabolomics", settings$id_column, settings$label_token)
  )

  multiomics <- collect_multiomics_layers(layers)
  multiomics$config <- config
  multiomics
}
