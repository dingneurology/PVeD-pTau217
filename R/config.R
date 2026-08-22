load_analysis_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install it with install.packages('yaml').", call. = FALSE)
  }
  config_path <- normalizePath(path, mustWork = TRUE)
  config <- yaml::read_yaml(config_path)
  if (!is.list(config)) {
    stop("Configuration root must be a mapping.", call. = FALSE)
  }
  list(
    config = config,
    project_root = normalizePath(file.path(dirname(config_path), ".."), mustWork = TRUE)
  )
}

resolve_project_path <- function(project_root, configured_path) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", configured_path)) {
    normalizePath(configured_path, mustWork = FALSE)
  } else {
    normalizePath(file.path(project_root, configured_path), mustWork = FALSE)
  }
}

configured_path <- function(config, project_root, key) {
  if (is.null(config$paths[[key]])) {
    stop(sprintf("Missing configuration path: paths.%s", key), call. = FALSE)
  }
  resolve_project_path(project_root, config$paths[[key]])
}
