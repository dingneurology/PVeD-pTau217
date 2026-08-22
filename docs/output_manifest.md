# Manuscript output manifest

This manifest prevents a repository from being described as fully reproducible
before every reported item has a single generating path.

| Manuscript item | Generating code | Public input | Current status |
|---|---|---|---|
| Cohort matching counts | `02_build_habs_time_aligned.R` | none; aggregate output after review | Implemented |
| ADNI primary interaction table | `03_run_primary_models.R` | controlled local table | Implemented; numeric rerun required |
| HABS primary/sensitivity interaction table | `02` then `03` | controlled local tables | Implemented; date verification and rerun required |
| ADNI and HABS longitudinal cognition table | `04_run_habs_cognition.R` | controlled local tables | Implemented; cognitive-input validation required |
| Cross-cohort interaction forest | `05_make_figures.R` | aggregate model results | Implemented |
| Manuscript Figure 1 cohort/profile layout | not yet assigned | disclosure-reviewed aggregate counts | Pending final manuscript layout |
| Manuscript Figure 2 phenotype panels | not yet assigned | disclosure-reviewed aggregate summaries | Pending |
| Manuscript Figure 3 interaction evidence | `05_make_figures.R` for forest component | aggregate model results | Partly implemented |
| Manuscript Figure 4 longitudinal cognition visual | not yet assigned | aggregate mixed-model contrasts | Pending final model |
| Manuscript Figure 4 longitudinal cognition visual | aggregate profile-slope results | aggregate model results | Source table implemented; final plotting script pending |
| Supplementary tables and figures | not yet frozen | aggregate outputs | Pending statistical-analysis-plan lock |

The pending rows are release blockers for a claim of complete
figure-to-code reproducibility. They are not reasons to include the old
participant-level source-data files in GitHub.
