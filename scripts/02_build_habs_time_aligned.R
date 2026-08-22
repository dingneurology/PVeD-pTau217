#!/usr/bin/env Rscript

# Builds an all-available primary HABS-HD table and a ±24-month sensitivity table.
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 2 && args[[1]] == "--config") args[[2]] else "config/analysis.yml"

source(file.path("R", "config.R"))
source(file.path("R", "validation.R"))
source(file.path("R", "time_alignment.R"))

loaded <- load_analysis_config(config_path)
config <- loaded$config
project_root <- loaded$project_root
participant <- config$columns$participant_id
pved <- config$columns$pved
pved_scan_date <- config$columns$pved_scan_date
measurement_date <- config$columns$measurement_date
tie_break <- config$time_matching$tie_break

index <- utils::read.csv(configured_path(config, project_root, "habs_index"), stringsAsFactors = FALSE)
covariates <- utils::read.csv(configured_path(config, project_root, "habs_covariates"), stringsAsFactors = FALSE)
ptau <- utils::read.csv(configured_path(config, project_root, "habs_ptau"), stringsAsFactors = FALSE)
amyloid <- utils::read.csv(configured_path(config, project_root, "habs_amyloid_pet"), stringsAsFactors = FALSE)
tau <- utils::read.csv(configured_path(config, project_root, "habs_tau_pet"), stringsAsFactors = FALSE)

require_columns(index, c(participant, pved, pved_scan_date), "HABS PVeD index")
require_unique(index, participant, "HABS PVeD index")
require_unique(covariates, participant, "HABS covariates")
index <- parse_required_dates(index, pved_scan_date, "HABS PVeD index", allow_missing = FALSE)

outcomes <- config$models$habs_outcomes
amyloid_outcomes <- names(outcomes)[vapply(outcomes, function(x) identical(x$family, "amyloid"), logical(1))]
tau_outcomes <- names(outcomes)[vapply(outcomes, function(x) identical(x$family, "tau"), logical(1))]
if (length(amyloid_outcomes) == 0 || length(tau_outcomes) == 0) stop("HABS outcome configuration requires amyloid and tau families.", call. = FALSE)

rename_if_present <- function(data, old, new) {
  if (old %in% names(data)) names(data)[names(data) == old] <- new
  data
}

build_habs_table <- function(window_months) {
  matched <- list(
    ptau_primary = closest_row_within_window(index, ptau, participant, pved_scan_date, measurement_date, config$columns$habs_ptau_primary, window_months, "ptau_primary", tie_break),
    ptau_alternative = closest_row_within_window(index, ptau, participant, pved_scan_date, measurement_date, config$columns$habs_ptau_alternative, window_months, "ptau_alternative", tie_break)
  )
  for (outcome in amyloid_outcomes) {
    matched[[paste0("amyloid_", outcome)]] <- closest_row_within_window(index, amyloid[, c(participant, measurement_date, outcome), drop = FALSE], participant, pved_scan_date, measurement_date, outcome, window_months, paste0("amyloid_", outcome), tie_break)
  }
  for (outcome in tau_outcomes) {
    matched[[paste0("tau_", outcome)]] <- closest_row_within_window(index, tau[, c(participant, measurement_date, outcome), drop = FALSE], participant, pved_scan_date, measurement_date, outcome, window_months, paste0("tau_", outcome), tie_break)
  }
  aligned <- Reduce(function(left, right) merge(left, right, by = participant, all.x = TRUE, sort = FALSE), c(list(index, covariates), matched))
  aligned <- rename_if_present(aligned, paste0("ptau_primary_", config$columns$habs_ptau_primary), config$columns$habs_ptau_primary)
  aligned <- rename_if_present(aligned, paste0("ptau_alternative_", config$columns$habs_ptau_alternative), config$columns$habs_ptau_alternative)
  for (outcome in names(outcomes)) {
    prefix <- if (identical(outcomes[[outcome]]$family, "amyloid")) "amyloid" else "tau"
    aligned <- rename_if_present(aligned, paste(prefix, outcome, outcome, sep = "_"), outcome)
  }
  aligned
}

all_available <- build_habs_table(NULL)
time_aligned <- build_habs_table(as.numeric(config$time_matching$sensitivity_window_months))
write_csv_public(all_available, configured_path(config, project_root, "habs_all_available"))
write_csv_public(time_aligned, configured_path(config, project_root, "habs_time_aligned"))

audit_counts <- function(data, analysis_sample) {
  data.frame(
    analysis_sample = analysis_sample,
    stage = c("verified_pved_index", "primary_ptau_available", "alternative_ptau_available", "amyloid_pet_available", "tau_pet_available"),
    n = c(nrow(data), sum(!is.na(data[[config$columns$habs_ptau_primary]])), sum(!is.na(data[[config$columns$habs_ptau_alternative]])), sum(rowSums(!is.na(data[, amyloid_outcomes, drop = FALSE])) > 0), sum(rowSums(!is.na(data[, tau_outcomes, drop = FALSE])) > 0))
  )
}
write_csv_public(rbind(audit_counts(all_available, "all_available"), audit_counts(time_aligned, "time_aligned_24m")), file.path(configured_path(config, project_root, "aggregate_results"), "habs_matching_counts.csv"))
message("WROTE HABS ALL-AVAILABLE PRIMARY AND ±24-MONTH SENSITIVITY TABLES")
message("No participant-level file was written to a public results directory.")
