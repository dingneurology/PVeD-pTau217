source(file.path("R", "validation.R"))

run_profile_slope_model <- function(baseline, cognition, cohort, participant, visit_date, outcome, ptau, pved, covariates, categorical, minimum_visits = 2) {
  require_columns(baseline, c(participant, ptau, pved, covariates), paste(cohort, "baseline"))
  require_columns(cognition, c(participant, visit_date, outcome), paste(cohort, "cognition"))
  visits <- merge(cognition, baseline[, c(participant, ptau, pved, covariates), drop = FALSE], by = participant, all = FALSE)
  visits <- parse_required_dates(visits, visit_date, paste(cohort, "cognition"), allow_missing = FALSE)
  visits <- visits[stats::complete.cases(visits[, c(outcome, ptau, pved, covariates), drop = FALSE]), , drop = FALSE]
  visits[[outcome]] <- as.numeric(visits[[outcome]])
  visits$time_years <- ave(as.numeric(visits[[visit_date]]), visits[[participant]], FUN = function(x) (x - min(x)) / 365.25)
  counts <- table(visits[[participant]])
  visits <- visits[visits[[participant]] %in% names(counts[counts >= minimum_visits]), , drop = FALSE]
  slopes <- do.call(rbind, lapply(split(visits, visits[[participant]]), function(x) {
    fit <- stats::lm(stats::as.formula(paste(outcome, "~ time_years")), data = x)
    data.frame(slope = unname(stats::coef(fit)[["time_years"]]), x[1, c(participant, ptau, pved, covariates), drop = FALSE])
  }))
  slopes$profile <- paste(ifelse(slopes[[ptau]] >= stats::median(slopes[[ptau]]), "high_pTau217", "low_pTau217"), ifelse(slopes[[pved]] >= stats::median(slopes[[pved]]), "high_PVeD", "low_PVeD"), sep = "_")
  slopes$profile <- stats::relevel(as.factor(slopes$profile), ref = "high_pTau217_high_PVeD")
  for (covariate in intersect(covariates, categorical)) slopes[[covariate]] <- as.factor(slopes[[covariate]])
  fit <- stats::lm(stats::as.formula(paste("slope ~", paste(c("profile", covariates), collapse = " + "))), data = slopes)
  coefficient <- summary(fit)$coefficients["profilehigh_pTau217_low_PVeD", , drop = FALSE]
  data.frame(cohort = cohort, model_id = "profile_slope_primary", contrast = "high_pTau217_low_PVeD_vs_high_pTau217_high_PVeD", n_participants = nrow(slopes), estimate = coefficient[1, "Estimate"], standard_error = coefficient[1, "Std. Error"], p_value = coefficient[1, "Pr(>|t|)"], covariates = paste(covariates, collapse = ";"))
}
