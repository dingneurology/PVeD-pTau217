# Reproducibility

## 1. Environment

- Primary language: R 4.4.1 (the version reported in the manuscript)
- R dependencies: `yaml`, `ggplot2`, and `lme4`
- Auxiliary language: Python 3.12.4 for the public release audit
- Random seed: `217` unless an analysis explicitly declares another seed
- Dates: ISO 8601 (`YYYY-MM-DD`)
- Time interval: days divided by 365.25 for years and by 30.4375 for months

Create a clean environment and install dependencies:

```bash
Rscript -e 'install.packages(c("yaml", "ggplot2", "lme4"), repos = "https://cloud.r-project.org")'
Rscript tests/r/run_tests.R

python scripts/00_release_audit.py
```

## 2. Controlled inputs

Copy `config/analysis.example.yml` to the ignored file
`config/analysis.yml`. Place approved source tables under `data/raw/`, or point
the configuration to approved local paths outside the repository. Never commit
the local configuration when it contains absolute paths.

Input contracts are documented in `docs/data_dictionary.md`. The pipeline
requires an exact PVeD scan date for HABS-HD. It intentionally stops if that
date is missing; it does not infer the date from the earliest available MRI.

## 3. Ordered execution

```bash
python scripts/00_release_audit.py
Rscript scripts/01_validate_inputs.R --config config/analysis.yml
Rscript scripts/02_build_habs_time_aligned.R --config config/analysis.yml
Rscript scripts/03_run_primary_models.R --config config/analysis.yml
Rscript scripts/04_run_habs_cognition.R --config config/analysis.yml
Rscript scripts/05_make_figures.R --config config/analysis.yml
python scripts/00_release_audit.py
```

The HABS time-aligned table is participant-level and is written only to the
ignored `data/derived/` directory.

## 4. Expected public outputs

Before manuscript submission, disclosure-reviewed aggregate outputs may be
copied deliberately into `results/public/`:

- model coefficients, robust standard errors, confidence intervals, P values,
  and FDR-adjusted P values;
- cohort flow counts that satisfy disclosure rules;
- aggregate plot source data without identifiers or exact dates;
- software/session information and an output manifest.

Generated participant-level tables must not be moved into a public directory.

## 5. Verification gates

A release is acceptable only when:

1. all synthetic tests pass;
2. all controlled input schemas pass;
3. HABS date matching uses verified PVeD scan dates and closest eligible
   measurements;
4. all model matrices are full rank and within the configured condition-number
   limit;
5. manuscript estimates match regenerated aggregate outputs;
6. `scripts/00_release_audit.py` exits successfully;
7. a human reviewer confirms the archive contains no controlled data;
8. the archive hash and Git commit are recorded in the manuscript package.

## 6. Known boundary

The repository cannot reproduce numeric results without separately authorised
cohort data. This is a controlled-data constraint, not permission to publish
model-ready participant-level files.
