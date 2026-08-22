# Public release checklist

## Scientific

- [ ] Primary HABS pTau217 assay is confirmed and matches the manuscript.
- [ ] Exact PVeD scan dates are verified.
- [ ] Outcome families and FDR scope are frozen.
- [ ] Models cannot include the outcome as a covariate.
- [ ] Longitudinal visits are anchored to the index date.
- [ ] Regenerated estimates match manuscript tables and figures.
- [ ] Claims remain associational and respect cohort/assay differences.

## Privacy and governance

- [ ] `data/raw/` and `data/derived/` contain no tracked files.
- [ ] No RID, PTID, Med_ID, image ID, exact visit date, or individual slope is
      present in public files.
- [ ] Aggregate source data satisfy cohort-specific small-cell rules.
- [ ] Data-use agreements and institutional requirements have been reviewed.
- [ ] Git history has been checked, not only the current working tree.

## Reproducibility

- [ ] Fresh-environment installation succeeds.
- [ ] Synthetic tests pass.
- [ ] Input validation passes on controlled data.
- [ ] All outputs have a generating script.
- [ ] Dependency versions and random seeds are recorded.
- [ ] Release audit passes.
- [ ] CITATION placeholders are replaced.
- [ ] Repository URL, article title, authors, and journal are final.
- [ ] Release commit and archive SHA-256 are recorded.

## Human sign-off

- [ ] Analysis lead
- [ ] Data governance lead
- [ ] Corresponding author
