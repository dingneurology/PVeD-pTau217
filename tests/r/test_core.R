set.seed(217)

source(file.path("R", "validation.R"))
source(file.path("R", "time_alignment.R"))
source(file.path("R", "models.R"))

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

index <- data.frame(
  participant_id = c("SYN001", "SYN002"),
  pved_scan_date = c("2020-01-15", NA),
  stringsAsFactors = FALSE
)
measurements <- data.frame(
  participant_id = c("SYN001", "SYN001"),
  measurement_date = c("2018-02-01", "2020-01-10"),
  ptau217_primary = c(1.1, 2.2),
  stringsAsFactors = FALSE
)

date_error <- tryCatch(
  closest_row_within_window(index, measurements, "participant_id", "pved_scan_date", "measurement_date", "ptau217_primary", 24, "ptau"),
  error = function(e) conditionMessage(e)
)
assert(grepl("missing or invalid dates", date_error), "Missing PVeD scan dates must fail.")

matched <- closest_row_within_window(
  index[1, , drop = FALSE],
  measurements,
  "participant_id",
  "pved_scan_date",
  "measurement_date",
  "ptau217_primary",
  24,
  "ptau"
)
assert(abs(matched$ptau_ptau217_primary - 2.2) < 1e-10, "Matching must select the closest eligible date, not the earliest date.")

n <- 120
pved <- rnorm(n)
ptau <- exp(rnorm(n))
ptau_z <- as.numeric(scale(log10(ptau)))
pved_z <- as.numeric(scale(pved))
interaction <- ptau_z * pved_z
model_data <- data.frame(
  participant_id = sprintf("SYN%03d", seq_len(n)),
  ptau217 = ptau,
  pved = pved,
  age = rnorm(n, 72, 6),
  sex = sample(c("female", "male"), n, replace = TRUE),
  education = rnorm(n, 16, 2),
  amyloid_centiloid = -0.45 * interaction + rnorm(n, sd = 0.65),
  stringsAsFactors = FALSE
)

fit <- fit_interaction_ols(
  model_data,
  outcome = "amyloid_centiloid",
  ptau = "ptau217",
  pved = "pved",
  covariates = c("age", "sex", "education"),
  categorical = c("sex"),
  minimum_n = 30
)
assert(fit$n == n, "Complete-case sample size is incorrect.")
assert(fit$interaction_beta < 0, "Synthetic interaction direction was not recovered.")

bad_model <- tryCatch(
  fit_interaction_ols(
    model_data,
    outcome = "amyloid_centiloid",
    ptau = "ptau217",
    pved = "pved",
    covariates = c("age", "amyloid_centiloid"),
    categorical = character(),
    minimum_n = 30
  ),
  error = function(e) conditionMessage(e)
)
assert(grepl("also listed as a covariate", bad_model), "Outcome-as-covariate models must be rejected.")

cat("R CORE TESTS PASSED\n")
