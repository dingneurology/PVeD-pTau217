# Legacy-to-release code map

The original project is preserved outside this public repository. It contains
exploratory scripts, private absolute paths, duplicate generations of some
tables and figures, and participant-level exports. Copying those files directly
would preserve known defects and create a disclosure risk.

The release candidate therefore uses a clean, reviewed implementation:

| Original analytical purpose | Public release entry point | Resolution |
|---|---|---|
| ADNI continuous interaction tables | `scripts/03_run_primary_models.R` | R `lm()` with HC3 inference, explicit complete cases, rank checks |
| HABS time-window validation | `scripts/02_build_habs_time_aligned.R` | verified index date and closest eligible measurement |
| HABS primary interaction tables | `scripts/03_run_primary_models.R` | all-available primary analysis, temporal and assay sensitivities |
| HABS PACC-like longitudinal validation | `scripts/04_run_habs_cognition.R` | R mixed-effects model; post-index visits; no individual-slope export |
| Main interaction forest figure | `scripts/05_make_figures.R` | reads aggregate model output only |
| Public-package validation | `scripts/00_release_audit.py` | auxiliary Python audit blocks controlled files and identifier-bearing tables |

## Deliberately retired behaviours

- hard-coded paths under an individual user's home directory;
- assignment of the earliest MRI as the PVeD scan date;
- selection of the earliest pTau/PET observation before time-window filtering;
- categorical encoding before missing-data filtering;
- using the outcome itself as an adjustment covariate;
- manual normal-approximation inference without design diagnostics;
- ambiguous primary-versus-alternative assay labels;
- public export of participant-level cognitive visits or slopes.

The previous Python implementation is retained under `python/` as auxiliary
reference material. The manuscript-facing analysis path is R-first. This map
documents analytical continuity, but numeric equivalence is neither assumed nor
desired where the prior behaviour was defective. All revised results must be
regenerated and reconciled with the manuscript.
