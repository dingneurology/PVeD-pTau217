# Data availability

This repository contains analysis code, documentation, and synthetic tests. It
does not contain ADNI or HABS-HD participant-level data.

ADNI data must be requested from the Alzheimer's Disease Neuroimaging
Initiative under its current data-use terms. HABS-HD data must be requested from
the study investigators or designated repository under the applicable
data-use agreement.

The following must remain local and are excluded by `.gitignore`:

- raw cohort tables and imaging manifests;
- participant identifiers, including RID, PTID, Med_ID, image IDs, and scan IDs;
- dates and visit-level records;
- participant-level model-ready or time-aligned datasets;
- individual cognitive trajectories or individual slope estimates;
- intermediate files from joins, exclusions, or quality control.

Only disclosure-reviewed aggregate tables and source data may be considered for
public release. An aggregate file must:

1. contain no direct or study-specific participant identifier;
2. contain no participant-level rows or exact visit dates;
3. avoid small cells according to the governing data-use agreement;
4. be traceable to a script and a prespecified analysis;
5. pass `python scripts/00_release_audit.py`.

Passing the automated audit does not replace investigator review or compliance
with ADNI, HABS-HD, institutional, journal, or funder requirements.
