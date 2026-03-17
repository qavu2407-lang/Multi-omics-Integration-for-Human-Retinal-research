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

  sample_ids <- setdiff(names(feature_tbl), id_column)
  feature_ids <- feature_tbl[[id_column]]

  if (anyDuplicated(feature_ids)) {
    stop(sprintf("Duplicated feature identifiers detected in %s.", path))
  }

  mat <- feature_tbl %>%
    tibble::column_to_rownames(id_column) %>%
    as.matrix()
  storage.mode(mat) <- "double"

  sample_metadata <- tibble::tibble(
    sample_id = sample_ids,
    group = as.character(unlist(label_row[1, sample_ids, drop = TRUE])),
    view = view_name
  )

  list(
    view = view_name,
    matrix = mat,
    sample_metadata = sample_metadata,
    raw_dimensions = dim(raw_tbl)
  )
}

align_multiomics_layers <- function(layer_list) {
  sample_sets <- purrr::map(layer_list, ~ colnames(.x$matrix))
  common_samples <- Reduce(intersect, sample_sets)

  if (length(common_samples) == 0) {
    stop("No shared samples were found across omics layers.")
  }

  aligned_layers <- purrr::map(layer_list, function(layer) {
    idx <- match(common_samples, colnames(layer$matrix))
    layer$matrix <- layer$matrix[, idx, drop = FALSE]
    layer$sample_metadata <- layer$sample_metadata %>%
      dplyr::filter(sample_id %in% common_samples) %>%
      dplyr::mutate(sample_id = factor(sample_id, levels = common_samples)) %>%
      dplyr::arrange(sample_id) %>%
      dplyr::mutate(sample_id = as.character(sample_id))
    layer
  })

  reference_groups <- aligned_layers[[1]]$sample_metadata %>%
    dplyr::select(sample_id, group)

  group_consistency <- purrr::map_lgl(aligned_layers[-1], function(layer) {
    current <- layer$sample_metadata %>%
      dplyr::select(sample_id, group)
    identical(reference_groups, current)
  })

  if (!all(group_consistency)) {
    stop("Sample group labels are not consistent across omics layers.")
  }

  list(
    views = purrr::set_names(purrr::map(aligned_layers, "matrix"), purrr::map_chr(aligned_layers, "view")),
    sample_metadata = reference_groups
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

  aligned <- align_multiomics_layers(layers)
  aligned$config <- config
  aligned
}
