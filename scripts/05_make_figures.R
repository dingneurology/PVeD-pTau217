#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 2 && args[[1]] == "--config") args[[2]] else "config/analysis.yml"

source(file.path("R", "config.R"))
source(file.path("R", "figures.R"))

loaded <- load_analysis_config(config_path)
config <- loaded$config
project_root <- loaded$project_root

results_path <- file.path(configured_path(config, project_root, "aggregate_results"), "primary_interaction_models.csv")
results <- utils::read.csv(results_path, stringsAsFactors = FALSE)
output_path <- file.path(project_root, "figures", "interaction_forest.png")
make_interaction_forest(results, output_path)
message(sprintf("WROTE FIGURE: %s", output_path))
