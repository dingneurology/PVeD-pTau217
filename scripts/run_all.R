#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 2 && args[[1]] == "--config") args[[2]] else "config/analysis.yml"

steps <- c(
  "scripts/01_validate_inputs.R",
  "scripts/02_build_habs_time_aligned.R",
  "scripts/03_run_primary_models.R",
  "scripts/04_run_habs_cognition.R",
  "scripts/05_make_figures.R"
)

for (step in steps) {
  status <- system2("Rscript", c(step, "--config", config_path))
  if (status != 0) {
    stop(sprintf("Pipeline step failed: %s", step), call. = FALSE)
  }
}

status <- system2("python3", c("scripts/00_release_audit.py"))
if (status != 0) {
  stop("Release audit failed.", call. = FALSE)
}
