source(file.path("R", "validation.R"))

z_score <- function(x) {
  values <- suppressWarnings(as.numeric(x))
  sd_value <- stats::sd(values, na.rm = TRUE)
  if (!is.finite(sd_value) || sd_value == 0) {
    return(rep(NA_real_, length(values)))
  }
  (values - mean(values, na.rm = TRUE)) / sd_value
}

bh_fdr <- function(p_values) {
  out <- rep(NA_real_, length(p_values))
  ok <- !is.na(p_values) & is.finite(p_values)
  out[ok] <- stats::p.adjust(p_values[ok], method = "BH")
  out
}

bind_rows_fill <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  filled <- lapply(rows, function(row) {
    missing <- setdiff(all_names, names(row))
    for (column in missing) row[[column]] <- NA
    row[, all_names, drop = FALSE]
  })
  do.call(rbind, filled)
}

hc3_covariance <- function(fitted) {
  design <- stats::model.matrix(fitted)
  residuals <- stats::residuals(fitted)
  hat_values <- stats::hatvalues(fitted)
  bread <- solve(crossprod(design))
  meat <- crossprod(design, design * as.vector((residuals / (1 - hat_values))^2))
  bread %*% meat %*% bread
}

prepare_interaction_design <- function(data,
                                       outcome,
                                       ptau,
                                       pved,
                                       covariates,
                                       categorical) {
  if (outcome %in% covariates) {
    stop(sprintf("Invalid model: outcome '%s' is also listed as a covariate.", outcome), call. = FALSE)
  }

  required <- c(outcome, ptau, pved, covariates)
  require_columns(data, required, "model table")
  require_positive(data, ptau, "model table", allow_missing = TRUE)

  raw <- data[, required, drop = FALSE]
  missing_by_variable <- vapply(raw, function(x) sum(is.na(x)), integer(1))
  complete <- raw[stats::complete.cases(raw), , drop = FALSE]
  complete[[ptau]] <- log10(as.numeric(complete[[ptau]]))

  numeric_for_check <- c(outcome, ptau, pved, setdiff(covariates, categorical))
  require_finite_variation(complete, numeric_for_check)

  design <- data.frame(
    outcome_z = z_score(complete[[outcome]]),
    ptau_z = z_score(complete[[ptau]]),
    pved_z = z_score(complete[[pved]])
  )
  design$ptau_x_pved <- design$ptau_z * design$pved_z

  for (covariate in covariates) {
    if (covariate %in% categorical) {
      factor_value <- droplevels(as.factor(complete[[covariate]]))
      if (nlevels(factor_value) < 2) {
        stop(sprintf("Categorical covariate '%s' has fewer than two levels.", covariate), call. = FALSE)
      }
      contrasts_matrix <- stats::model.matrix(~ factor_value)[, -1, drop = FALSE]
      colnames(contrasts_matrix) <- paste(covariate, levels(factor_value)[-1], sep = "_")
      design <- cbind(design, contrasts_matrix)
    } else {
      design[[covariate]] <- z_score(complete[[covariate]])
    }
  }

  if (anyNA(design)) {
    stop("Unexpected missing value after complete-case preparation.", call. = FALSE)
  }

  list(design = design, missing_by_variable = missing_by_variable)
}

fit_interaction_ols <- function(data,
                                outcome,
                                ptau,
                                pved,
                                covariates,
                                categorical,
                                minimum_n = 30,
                                max_condition_number = 1000000) {
  prepared <- prepare_interaction_design(data, outcome, ptau, pved, covariates, categorical)
  design <- prepared$design
  n <- nrow(design)
  if (n < minimum_n) {
    stop(sprintf("Complete-case sample is %d, below minimum %d.", n, minimum_n), call. = FALSE)
  }

  predictors <- setdiff(names(design), "outcome_z")
  model_formula <- stats::as.formula(paste("outcome_z ~", paste(predictors, collapse = " + ")))
  model_matrix <- stats::model.matrix(model_formula, data = design)
  matrix_rank <- qr(model_matrix)$rank
  if (matrix_rank != ncol(model_matrix)) {
    stop(sprintf("Design matrix is rank deficient: rank=%d, columns=%d.", matrix_rank, ncol(model_matrix)), call. = FALSE)
  }
  condition_number <- kappa(model_matrix, exact = TRUE)
  if (!is.finite(condition_number) || condition_number > max_condition_number) {
    stop(
      sprintf("Design condition number %.3g exceeds %.3g.", condition_number, max_condition_number),
      call. = FALSE
    )
  }

  fitted <- stats::lm(model_formula, data = design)
  covariance <- hc3_covariance(fitted)
  term <- "ptau_x_pved"
  beta <- stats::coef(fitted)[[term]]
  se <- sqrt(diag(covariance))[[term]]
  statistic <- beta / se
  p_value <- 2 * stats::pt(abs(statistic), df = stats::df.residual(fitted), lower.tail = FALSE)
  ci <- beta + c(-1, 1) * stats::qt(0.975, stats::df.residual(fitted)) * se

  list(
    n = n,
    interaction_beta = unname(beta),
    robust_se = unname(se),
    ci_low = unname(ci[[1]]),
    ci_high = unname(ci[[2]]),
    p_value = unname(p_value),
    r_squared = summary(fitted)$r.squared,
    adjusted_r_squared = summary(fitted)$adj.r.squared,
    covariance = "HC3",
    matrix_rank = matrix_rank,
    condition_number = condition_number,
    missing_by_variable = prepared$missing_by_variable
  )
}

run_outcome_family <- function(data,
                               cohort,
                               marker_role,
                               analysis_sample = "primary",
                               ptau,
                               pved,
                               outcomes,
                               model_sets,
                               categorical,
                               minimum_n,
                               max_condition_number) {
  rows <- list()
  row_index <- 1L
  for (outcome in names(outcomes)) {
    if (!outcome %in% names(data)) next
    metadata <- outcomes[[outcome]]
    for (model_id in names(model_sets)) {
      covariates <- unlist(model_sets[[model_id]], use.names = FALSE)
      base <- data.frame(
        cohort = cohort,
      marker_role = marker_role,
      analysis_sample = analysis_sample,
        outcome = outcome,
        outcome_label = metadata$label,
        analysis_family = metadata$family,
        model_id = model_id,
        covariates = paste(covariates, collapse = ";"),
        stringsAsFactors = FALSE
      )
      if (outcome %in% covariates) {
        rows[[row_index]] <- cbind(base, model_status = "skipped_outcome_is_covariate")
        row_index <- row_index + 1L
        next
      }
      result <- tryCatch(
        fit_interaction_ols(data, outcome, ptau, pved, covariates, categorical, minimum_n, max_condition_number),
        error = function(e) e
      )
      if (inherits(result, "error")) {
        rows[[row_index]] <- cbind(base, model_status = "failed_validation", note = conditionMessage(result))
      } else {
        rows[[row_index]] <- cbind(
          base,
          model_status = "ok",
          as.data.frame(result[setdiff(names(result), "missing_by_variable")], stringsAsFactors = FALSE)
        )
      }
      row_index <- row_index + 1L
    }
  }

  if (length(rows) == 0) return(data.frame())
  results <- bind_rows_fill(rows)
  results$fdr <- NA_real_
  ok <- results$model_status == "ok"
  if (any(ok)) {
    group_key <- paste(results$cohort, results$marker_role, results$analysis_sample, results$model_id, results$analysis_family, sep = "\r")
    for (key in unique(group_key[ok])) {
      subset_index <- ok & group_key == key
      results$fdr[subset_index] <- bh_fdr(as.numeric(results$p_value[subset_index]))
    }
  }
  results
}
