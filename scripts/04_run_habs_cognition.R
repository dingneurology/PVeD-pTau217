#!/usr/bin/env Rscript

# Manuscript-facing profile-slope contrasts for ADNI and HABS-HD.
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 2 && args[[1]] == "--config") args[[2]] else "config/analysis.yml"
source(file.path("R", "config.R"))
source(file.path("R", "validation.R"))
source(file.path("R", "cognition.R"))
loaded <- load_analysis_config(config_path)
config <- loaded$config
project_root <- loaded$project_root
participant <- config$columns$participant_id
categorical <- unlist(config$models$categorical_covariates, use.names = FALSE)
adni <- utils::read.csv(configured_path(config, project_root, "adni_model_ready"), stringsAsFactors = FALSE)
habs <- utils::read.csv(configured_path(config, project_root, "habs_all_available"), stringsAsFactors = FALSE)
adni_cognition <- utils::read.csv(configured_path(config, project_root, "adni_cognition"), stringsAsFactors = FALSE)
habs_cognition <- utils::read.csv(configured_path(config, project_root, "habs_cognition"), stringsAsFactors = FALSE)

adni_result <- run_profile_slope_model(adni, adni_cognition, "ADNI", participant, config$columns$cognition_visit_date, config$columns$adni_cognition_outcome, config$columns$adni_ptau, config$columns$pved, config$cognition$primary_covariates, categorical, config$cognition$minimum_visits)
habs_result <- run_profile_slope_model(habs, habs_cognition, "HABS-HD", participant, config$columns$cognition_visit_date, config$columns$cognition_outcome, config$columns$habs_ptau_primary, config$columns$pved, config$cognition$primary_covariates, categorical, config$cognition$minimum_visits)
output_dir <- configured_path(config, project_root, "aggregate_results")
write_csv_public(rbind(adni_result, habs_result), file.path(output_dir, "longitudinal_profile_slope_contrasts.csv"))
message("WROTE ADNI AND HABS-HD LONGITUDINAL PROFILE-SLOPE CONTRASTS")
