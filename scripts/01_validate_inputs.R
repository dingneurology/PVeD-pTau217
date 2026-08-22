#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 2 && args[[1]] == "--config") args[[2]] else "config/analysis.yml"

source(file.path("R", "config.R"))
source(file.path("R", "validation.R"))

loaded <- load_analysis_config(config_path)
config <- loaded$config
project_root <- loaded$project_root

required_paths <- c(
  "adni_model_ready",
  "habs_index",
  "habs_ptau",
  "habs_amyloid_pet",
  "habs_tau_pet",
  "habs_covariates",
  "habs_cognition",
  "adni_cognition"
)

missing_paths <- character()
for (key in required_paths) {
  path <- configured_path(config, project_root, key)
  if (!file.exists(path)) {
    missing_paths <- c(missing_paths, sprintf("paths.%s -> %s", key, path))
  }
}
if (length(missing_paths) > 0) {
  stop(sprintf("Required controlled input files are missing:\n%s", paste(missing_paths, collapse = "\n")), call. = FALSE)
}

participant <- config$columns$participant_id
pved <- config$columns$pved
pved_scan_date <- config$columns$pved_scan_date
measurement_date <- config$columns$measurement_date

habs_index <- utils::read.csv(configured_path(config, project_root, "habs_index"), stringsAsFactors = FALSE)
require_columns(habs_index, c(participant, pved, pved_scan_date), "HABS PVeD index")
require_unique(habs_index, participant, "HABS PVeD index")
parse_required_dates(habs_index, pved_scan_date, "HABS PVeD index", allow_missing = FALSE)

habs_ptau <- utils::read.csv(configured_path(config, project_root, "habs_ptau"), stringsAsFactors = FALSE)
require_columns(habs_ptau, c(participant, measurement_date, config$columns$habs_ptau_primary, config$columns$habs_ptau_alternative), "HABS pTau217")
parse_required_dates(habs_ptau, measurement_date, "HABS pTau217", allow_missing = TRUE)

adni <- utils::read.csv(configured_path(config, project_root, "adni_model_ready"), stringsAsFactors = FALSE)
require_columns(adni, c(participant, pved, config$columns$adni_ptau, config$models$adni$primary_covariates), "ADNI model-ready table")

message("INPUT VALIDATION PASSED")
message("Exact HABS PVeD scan dates are present; no earliest-MRI fallback is allowed.")
