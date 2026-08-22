#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 2 && args[[1]] == "--config") args[[2]] else "config/analysis.yml"

source(file.path("R", "config.R"))
source(file.path("R", "validation.R"))
source(file.path("R", "models.R"))

loaded <- load_analysis_config(config_path)
config <- loaded$config
project_root <- loaded$project_root
model_sets <- function(section) c(list(primary = section$primary_covariates), section$sensitivity_covariates)
categorical <- unlist(config$models$categorical_covariates, use.names = FALSE)
minimum_n <- as.integer(config$project$minimum_complete_n)
max_condition_number <- as.numeric(config$project$max_condition_number)

adni <- utils::read.csv(configured_path(config, project_root, "adni_model_ready"), stringsAsFactors = FALSE)
habs_all_available <- utils::read.csv(configured_path(config, project_root, "habs_all_available"), stringsAsFactors = FALSE)
habs_time_aligned <- utils::read.csv(configured_path(config, project_root, "habs_time_aligned"), stringsAsFactors = FALSE)
pved <- config$columns$pved

adni_results <- run_outcome_family(adni, cohort = "ADNI", marker_role = "discovery", analysis_sample = "primary", ptau = config$columns$adni_ptau, pved = pved, outcomes = config$models$adni_outcomes, model_sets = model_sets(config$models$adni), categorical = categorical, minimum_n = minimum_n, max_condition_number = max_condition_number)
habs_primary <- run_outcome_family(habs_all_available, cohort = "HABS-HD", marker_role = "primary", analysis_sample = "all_available", ptau = config$columns$habs_ptau_primary, pved = pved, outcomes = config$models$habs_outcomes, model_sets = model_sets(config$models$habs), categorical = categorical, minimum_n = minimum_n, max_condition_number = max_condition_number)
habs_time_sensitivity <- run_outcome_family(habs_time_aligned, cohort = "HABS-HD", marker_role = "primary", analysis_sample = "time_aligned_24m", ptau = config$columns$habs_ptau_primary, pved = pved, outcomes = config$models$habs_outcomes, model_sets = model_sets(config$models$habs), categorical = categorical, minimum_n = minimum_n, max_condition_number = max_condition_number)
habs_assay_sensitivity <- run_outcome_family(habs_all_available, cohort = "HABS-HD", marker_role = "alternative_assay", analysis_sample = "all_available", ptau = config$columns$habs_ptau_alternative, pved = pved, outcomes = config$models$habs_outcomes, model_sets = model_sets(config$models$habs), categorical = categorical, minimum_n = minimum_n, max_condition_number = max_condition_number)

results <- rbind(adni_results, habs_primary, habs_time_sensitivity, habs_assay_sensitivity)
output_dir <- configured_path(config, project_root, "aggregate_results")
write_csv_public(results, file.path(output_dir, "primary_interaction_models.csv"))
message(sprintf("WROTE AGGREGATE MODEL RESULTS: %s", file.path(output_dir, "primary_interaction_models.csv")))
