# Input data contract

All tables described here are controlled local inputs and must remain outside
Git. Column names below are canonical names used by the example configuration.

## ADNI model-ready table

One row per participant.

Required:

- `participant_id`: local study identifier; never exported;
- `ptau217`: positive plasma pTau217 concentration;
- `pved`: continuous PVeD measure;
- `age`, `sex`, `education`;
- one or more configured ADNI outcomes.

Optional sensitivity variables include `apoe4`.

## HABS PVeD index table

Exactly one row per participant:

- `participant_id`;
- `pved`;
- `pved_scan_date`: verified date of the MRI used to derive PVeD.

## HABS measurement tables

Long format, with zero or more rows per participant:

- pTau table: `participant_id`, `measurement_date`, and the two study-specific
  columns mapped in the private configuration to `habs_ptau_primary` and
  `habs_ptau_alternative`;
- amyloid PET table: `participant_id`, `measurement_date`,
  `amyloid_centiloid`, `amyloid_suvr`;
- tau PET table: `participant_id`, `measurement_date` and configured tau
  outcomes.

Each row must represent a single measurement session. Do not collapse to the
earliest visit before running the matching script.

## HABS covariates

Exactly one prespecified baseline row per participant:

- `participant_id`, `age`, `sex`, `education`;
- optional `apoe4`.

The baseline definition must be documented upstream. Selecting the first
nonmissing value independently for each covariate is not permitted.

## HABS cognition

Long format:

- `participant_id`;
- `visit_date`;
- `pacc_like`;
- component-availability indicators recommended for audit.

Construction and norming of PACC-like must be documented separately. Combining
scores standardised against incompatible reference populations requires
validation.

## Missing values and dates

Convert study-specific missing codes to true missing values before these
tables are used. Dates must be valid ISO 8601 values. The pipeline rejects
duplicate index rows, missing index dates, non-positive pTau217 values, and
missing required columns.
