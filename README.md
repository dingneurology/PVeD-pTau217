# PVeD–pTau217 analysis

Reproducible analysis code for evaluating whether MRI-derived periventricular
diffusivity (PVeD) modifies associations of plasma pTau217 with Alzheimer disease
PET biomarkers and longitudinal cognition in ADNI and HABS-HD.

> **Release status:** review candidate. Statistical results must be regenerated
> after the controlled-data inputs and the exact PVeD scan dates have been
> verified. No participant-level data are included.

## Scientific scope

The primary estimand is the coefficient of the continuous interaction
`z(log10[pTau217]) × z(PVeD)` in prespecified outcome models. ADNI is the
discovery cohort. HABS-HD is an external validation cohort. Profile plots based
on median splits are descriptive and are not substitutes for the continuous
interaction analysis.

This repository supports association analyses. It does not establish causality,
clinical utility, or a direct measure of glymphatic function.

## Repository structure

```text
R/                      Reusable R functions for validation, matching, models, figures
config/                 Example analysis configuration
data/raw/               Local controlled inputs; ignored by Git
data/derived/           Local participant-level derivatives; ignored by Git
docs/                   Analysis plan, data contract, and release checklist
figures/                Generated figures; ignored until publication approval
src/pved_ptau217/       Python helper used by the release privacy audit
results/                Generated aggregate results; ignored until approval
scripts/                R-first ordered command-line entry points
tests/                  Synthetic-data R unit tests
```

## Data boundary

ADNI and HABS-HD participant-level data are not redistributed. Investigators
must obtain data through the corresponding data-use procedures. Place approved
local inputs under `data/raw/`; Git ignores this directory.

Do not place names, dates of birth, medical-record identifiers, ADNI RID/PTID,
HABS Med_ID, image identifiers, or participant-level model-ready tables in a
public repository. See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).

## Quick start

R is the primary analysis language. Python is used only for the auxiliary public
release audit.

```bash
Rscript -e 'install.packages(c("yaml", "ggplot2", "lme4"), repos = "https://cloud.r-project.org")'
Rscript tests/r/run_tests.R

cp config/analysis.example.yml config/analysis.yml
python scripts/00_release_audit.py
Rscript scripts/01_validate_inputs.R --config config/analysis.yml
Rscript scripts/02_build_habs_time_aligned.R --config config/analysis.yml
Rscript scripts/03_run_primary_models.R --config config/analysis.yml
Rscript scripts/04_run_habs_cognition.R --config config/analysis.yml
Rscript scripts/05_make_figures.R --config config/analysis.yml
python scripts/00_release_audit.py
```

`config/analysis.yml` is intentionally ignored because local paths may disclose
private information. The example configuration uses repository-relative paths.

For a one-command local rerun after the configuration is complete:

```bash
Rscript scripts/run_all.R --config config/analysis.yml
```

## Analysis safeguards

- Complete cases are defined on the original variables before categorical
  encoding; missing categories cannot silently become reference levels.
- A model is rejected if the outcome is also listed as a covariate.
- HABS measurements are selected by minimum absolute distance from the verified
  PVeD scan date within the prespecified window.
- No fallback to the earliest MRI or earliest biomarker measurement is allowed.
- OLS inference uses R `lm()` with explicitly computed HC3 robust standard
  errors.
- Rank-deficient and severely ill-conditioned designs fail explicitly.
- False-discovery-rate families are declared in the configuration and retained
  in every output row.
- Longitudinal cognition uses aggregate participant-specific annual slopes for
  manuscript-facing profile contrasts; no participant-level slopes are exported.
- Publication figures are generated from aggregate model outputs rather than
  participant-level tables.
- The release audit blocks participant identifiers, participant-level tabular
  files, hard-coded home-directory paths, likely secrets, and controlled data.

## Reproduction

Detailed commands, expected inputs, and validation gates are in
[REPRODUCIBILITY.md](REPRODUCIBILITY.md). The prespecified estimands and
multiplicity families are in [docs/analysis_plan.md](docs/analysis_plan.md).

## Citation and licence

Citation metadata are provided in [CITATION.cff](CITATION.cff). The code is
released under the BSD 3-Clause License. Cohort data remain subject to their
original data-use agreements.
