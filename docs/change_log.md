# Change log

## 0.1.0 — 2026-07-27

- Rebuilt the public package with a code-only whitelist.
- Reorganised the manuscript-facing workflow as R-first code, with Python used
  only for auxiliary release auditing and legacy reference material.
- Removed all hard-coded user paths.
- Added explicit controlled-data boundaries.
- Added exact-date, closest-within-window HABS matching.
- Defined a primary HABS pTau217 assay and an alternative-assay sensitivity
  analysis without assay-version labels.
- Added complete-case filtering before categorical encoding.
- Added outcome-as-covariate rejection.
- Replaced manual normal-approximation OLS inference with R `lm()` and
  explicitly computed HC3 robust covariance.
- Added rank and condition-number checks.
- Replaced public individual-slope exports with aggregate mixed-model results.
- Added synthetic tests and a release privacy audit.

Numeric results have not yet been regenerated in this release candidate.
