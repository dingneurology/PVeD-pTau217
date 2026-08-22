# Code–manuscript alignment record

This release candidate is aligned to the submitted manuscript dated 2026-08-22.

- Terminology uses `pTau217`, `primary`, and `alternative_assay`; public files do not use PLUS or V2 assay labels.
- Baseline clinical diagnosis is not included as a covariate.
- ADNI primary cross-sectional models use age, sex, and education; HABS-HD primary cross-sectional models use age, sex, and APOE ε4.
- HABS-HD all-available analyses are primary; the ±24-month table is a temporal-alignment sensitivity analysis.
- The longitudinal manuscript-facing output is the adjusted annual-slope contrast between high-pTau217/low-PVeD and high-pTau217/high-PVeD profiles in both cohorts.
- This repository analyzes precomputed PVeD. Exact PVeD generation must be documented with the linked EstPVeD workflow and the study-specific preprocessing record.

Before pushing this release, run the full pipeline on controlled data, reconcile every reported estimate with the manuscript, and confirm the declared FDR families reproduce the submitted adjusted P values.
