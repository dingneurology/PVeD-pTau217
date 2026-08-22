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
