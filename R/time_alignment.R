closest_row_within_window <- function(index,
                                      measurements,
                                      participant,
                                      index_date,
                                      measurement_date,
                                      required_value,
                                      window_months,
                                      prefix,
                                      tie_break = "earlier") {
  if (!tie_break %in% c("earlier", "later")) {
    stop("tie_break must be 'earlier' or 'later'.", call. = FALSE)
  }
  if (!is.null(window_months) && (!is.numeric(window_months) || length(window_months) != 1 || window_months <= 0)) {
    stop("window_months must be NULL (all available observations) or one positive number.", call. = FALSE)
  }

  require_columns(index, c(participant, index_date), "index")
  require_columns(measurements, c(participant, measurement_date, required_value), "measurements")
  require_unique(index, participant, "index")

  left <- parse_required_dates(index, index_date, "index", allow_missing = FALSE)
  right <- parse_required_dates(measurements, measurement_date, "measurements", allow_missing = TRUE)
  right <- right[!is.na(right[[measurement_date]]) & !is.na(right[[required_value]]), , drop = FALSE]

  candidates <- merge(
    left[, c(participant, index_date), drop = FALSE],
    right,
    by = participant,
    all.x = TRUE,
    sort = FALSE
  )
  candidates$signed_months <- as.numeric(candidates[[measurement_date]] - candidates[[index_date]]) / 30.4375
  candidates$absolute_months <- abs(candidates$signed_months)
  eligible <- candidates[!is.na(candidates$absolute_months), , drop = FALSE]
  if (!is.null(window_months)) {
    eligible <- eligible[eligible$absolute_months <= window_months, , drop = FALSE]
  }

  if (nrow(eligible) == 0) {
    empty <- data.frame(character(), stringsAsFactors = FALSE)
    names(empty) <- participant
    return(empty)
  }

  eligible$.__tie_order <- if (tie_break == "earlier") eligible$signed_months else -eligible$signed_months
  eligible <- eligible[order(
    eligible[[participant]],
    eligible$absolute_months,
    eligible$.__tie_order,
    eligible[[measurement_date]]
  ), , drop = FALSE]
  selected <- eligible[!duplicated(eligible[[participant]]), , drop = FALSE]

  protected <- c(participant, index_date, ".__tie_order")
  for (column in setdiff(names(selected), protected)) {
    names(selected)[names(selected) == column] <- paste(prefix, column, sep = "_")
  }
  selected[[index_date]] <- NULL
  selected$.__tie_order <- NULL
  rownames(selected) <- NULL
  selected
}
