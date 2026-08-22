# Manuscript-aligned analysis plan

## Primary question and estimand

The primary estimand is the standardized coefficient of the continuous
`z(log10[pTau217]) × z(PVeD)` interaction in cross-sectional outcome models.
PVeD is an MRI-derived periventricular diffusion measure and is not interpreted
as a direct measure of glymphatic flow, a causal mechanism, or a replacement
for plasma pTau217 or PET.

## Cohort-specific primary models

- ADNI discovery models adjust for age, sex, and education.
- HABS-HD external-validation models adjust for age, sex, and APOE ε4 status.
- Baseline clinical diagnosis is not a covariate in any model.
- The HABS-HD primary pTau217 assay is labelled `primary`; the second released
  pTau217 assay is retained only as an assay-specific sensitivity analysis.

## HABS-HD timing

The all-available sample is the HABS-HD primary analysis. For every participant
with a verified PVeD scan date, the closest available pTau217 or PET measure is
selected without a time restriction. A separate ±24-month table is generated
for the prespecified temporal-alignment sensitivity analysis. The pipeline
stops if the exact PVeD scan date is unavailable and never substitutes the
earliest MRI.

## Multiplicity and reporting

Benjamini-Hochberg false-discovery-rate adjustment is performed within each
declared cohort, marker role, analysis sample, model specification, and outcome
family. These labels are retained in every aggregate result row. Final FDR
families must be reconciled with the submitted tables before release.

## Longitudinal cognition

Median-defined pTau217–PVeD profiles are descriptive. The manuscript-facing
longitudinal analysis compares participant-specific annual cognitive slopes,
with the high-pTau217/low-PVeD versus high-pTau217/high-PVeD contrast adjusted
for age, sex, and education. APOE ε4 adjustment and direct repeated-measures
models are sensitivity analyses. ADNI and HABS-HD cognition scripts must be
run only after the cohort-specific cognitive input contracts are verified.

## PVeD derivation boundary

This repository analyses precomputed PVeD values. PVeD derivation from diffusion
MRI is implemented separately in EstPVeD (https://github.com/ChangleChen/EstPVeD).
The submitted analysis configuration must record the exact EstPVeD release or
commit, DSI Studio build, and preprocessing parameters used for the study.
