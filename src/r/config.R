read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  cfg$project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
  cfg
}

project_path <- function(...) {
  file.path(normalizePath(".", winslash = "/", mustWork = TRUE), ...)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

ensure_parent_dir <- function(path) {
  ensure_dir(dirname(path))
}

write_session_metadata <- function(path) {
  ensure_parent_dir(path)
  session_info <- utils::sessionInfo()
  payload <- list(
    timestamp = as.character(Sys.time()),
    r_version = as.character(getRversion()),
    platform = R.version$platform,
    packages = list(
      basePkgs = session_info$basePkgs,
      otherPkgs = if (is.null(session_info$otherPkgs)) NULL else names(session_info$otherPkgs),
      loadedOnly = if (is.null(session_info$loadedOnly)) NULL else names(session_info$loadedOnly)
    )
  )
  jsonlite::write_json(payload, path = path, auto_unbox = TRUE, pretty = TRUE)
  path
}
