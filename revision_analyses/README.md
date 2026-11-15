# Revision analyses

This is the entry point for the revision analysis code referred to in the BMC Medicine resubmission.
The source package uploaded on 15 November 2026 organizes the implementation in the root `analysis/` and `figures/` directories.

- [Analysis scripts](../analysis/): cross-sectional models, continuous and profile longitudinal models, specificity, classification and tau-tracer analyses.
- [Run the complete Python pipeline](../run_all.py).
- [Input paths and environment variables](../config.py).
- [Python dependencies](../requirements.txt).
- [Figure scripts](../figures/) and [R figures](../figures/R/make_figures.R).
- [Input files, execution instructions and output-to-manuscript mapping](../README.md).

Run the pipeline from the repository root. Controlled cohort inputs must be obtained separately through ADNI and HABS-HD; no participant-level files are supplied here.

The earlier `R/`, `scripts/` and `src/` implementation remains in the repository for provenance. Use the root README and `run_all.py` for this resubmission workflow.

This upload records the supplied code. Syntax checks do not constitute an independent reproduction of the submitted scientific results.
