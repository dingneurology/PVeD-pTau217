require_columns <- function(data, columns, table_name) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(
      sprintf("%s is missing required columns: %s", table_name, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

require_unique <- function(data, column, table_name) {
  duplicates <- data[[column]][duplicated(data[[column]])]
  if (length(duplicates) > 0) {
    stop(sprintf("%s has duplicate participant rows in '%s'.", table_name, column), call. = FALSE)
  }
  invisible(TRUE)
}

parse_required_dates <- function(data, columns, table_name, allow_missing = FALSE) {
  out <- data
  for (column in columns) {
    parsed <- as.Date(out[[column]])
    if (!allow_missing && any(is.na(parsed))) {
      stop(sprintf("%s has missing or invalid dates in '%s'.", table_name, column), call. = FALSE)
    }
    bad <- !is.na(out[[column]]) & is.na(parsed)
    if (any(bad)) {
      stop(sprintf("%s has invalid dates in '%s'.", table_name, column), call. = FALSE)
    }
    out[[column]] <- parsed
  }
  out
}

require_positive <- function(data, column, table_name, allow_missing = TRUE) {
  values <- suppressWarnings(as.numeric(data[[column]]))
  checked <- if (allow_missing) values[!is.na(values)] else values
  if (any(is.na(checked)) || any(checked <= 0)) {
    stop(sprintf("%s requires positive finite values in '%s'.", table_name, column), call. = FALSE)
  }
  invisible(TRUE)
}

require_finite_variation <- function(data, columns, table_name = "model table") {
  for (column in columns) {
    values <- suppressWarnings(as.numeric(data[[column]]))
    if (any(!is.finite(values))) {
      stop(sprintf("%s has non-finite values in '%s'.", table_name, column), call. = FALSE)
    }
    if (stats::sd(values) == 0 || is.na(stats::sd(values))) {
      stop(sprintf("%s has no variation in '%s'.", table_name, column), call. = FALSE)
    }
  }
  invisible(TRUE)
}

write_csv_public <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data, path, row.names = FALSE, na = "")
}
