# Code update: 15 November 2026

Added the author-supplied `PVeD-pTau217_github.zip` resubmission code to the existing repository.

- New Python analysis and figure entry points, path configuration, dependencies and an R figure script.
- Root README describes the uploaded workflow; `revision_analyses/` provides the path referenced in the submitted supplement.
- Previous source files and commit history are retained.
- No controlled input datasets, model-ready participant files or generated analysis outputs are uploaded.
- Source ZIP SHA-256: `61da6c38156dacbc4a16931ba518ec51e49265751c9988f477ed03bc27d5c20d`.
- Validation scope: syntax checks only; numeric manuscript results were not independently regenerated as part of this upload.

## Historical release-candidate notes

The notes below describe the earlier implementation and its release status at that time.

# Release-candidate notes

This package is structurally ready for GitHub review but should not be tagged as
a final reproducible release until the following study-specific items are
completed:

1. replace the generic study-author entry in `CITATION.cff` and add the final
   repository URL and article citation;
2. verify the BSD 3-Clause licence choice with all code owners;
3. map controlled cohort fields to the canonical input contract;
4. recover and verify the exact HABS PVeD scan dates;
5. regenerate numeric results using the corrected pipeline;
6. reconcile every manuscript number, table, figure, and multiplicity family;
7. obtain data-governance and corresponding-author sign-off.

The prior scripts are not copied into this public package because they contain
hard-coded private paths, duplicate analysis versions, and participant-level
export behaviour. They remain preserved in the original project for provenance.

## R-first reorganisation

The manuscript-facing analysis has been reorganised as an R-first workflow:

- reusable R functions are in `R/`;
- ordered article-logic entry points are in `scripts/*.R`;
- Python is retained only for the public release audit and auxiliary legacy
  reference material under `python/`;
- synthetic safeguards are tested through `Rscript tests/r/run_tests.R`.
