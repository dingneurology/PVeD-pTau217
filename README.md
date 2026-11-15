# PVeD-pTau217

Analysis code for:

> **Periventricular diffusivity modifies the association of plasma pTau217 with amyloid and tau PET: a discovery–validation study in ADNI and HABS-HD.**
> Bao R, Zhang S, He H, *et al.* (manuscript under revision at *BMC Medicine*).

The repository reproduces every analysis reported in the revised manuscript: the cross-sectional pTau217 × PVeD interaction models (Tables 2, S4–S13), the longitudinal mixed-effects models of PACC change (Tables 3, S14–S17, S21), the two-cohort synthesis (Table S22) and Figures 1 (data panels), 2, 3, 4 and Supplementary Figures S1–S3.

Periventricular diffusivity (PVeD) itself is computed with the unmodified EstPVeD pipeline: <https://github.com/ChangleChen/EstPVeD>.

## Data access

No participant-level data are included, and none may be committed to this repository.

| Cohort | Role | Access |
|---|---|---|
| ADNI | Discovery (n = 602) | <https://adni.loni.usc.edu/>, after registration and approval of the ADNI Data Use Agreement |
| HABS-HD | External validation (n = 2313) | <https://healthandagingbrainstudy.org/>, after application and approval |

The ADSP Phenotype Harmonization Consortium (PHC) tau PET file used for the tracer sensitivity analysis is distributed through ADNI/LONI.

## Repository layout

```
config.py                 input/output paths (edit here, or set environment variables)
run_all.py                runs the Python pipeline in order
requirements.txt
analysis/
  statlib.py              OLS, logistic regression, DeLong test, NRI/IDI, Benjamini-Hochberg,
                          REML linear mixed model (random intercept and slope)
  build_adni.py           ADNI: PVeD baseline file + PET, MRI and CDR closest in time to the PVeD MRI
  build_habs.py           HABS-HD: PVeD, plasma, PET, MRI, clinical, genomic and DTI-protocol releases
  prep.py                 harmonised analytic datasets, within-cohort standardisation, profiles
  analyses.py             cross-sectional models M1-M4, simple slopes, temporal alignment,
                          amyloid PET classification, specificity, robustness, Table 1, LMMs
  lmm_profile.py          profile-by-time mixed models of PACC change
  incr.py                 incremental explained variance
  extra.py                additional sensitivity analyses (amyloid-positive subgroup, hemisphere,
                          alternative HABS-HD assay, PACC definitions)
  tau_tracer.py           ADNI tau PET tracer handling (flortaucipir-only primary analysis)
  long_methods.py         ADNI longitudinal model-form comparison
  meta_long.py            two-cohort DerSimonian-Laird synthesis
figures/
  figs.py                 Figure 4, Supplementary Figures S2-S3 (matplotlib)
  export_fig_sources.py   source data for the R figures
  R/make_figures.R        Figures 2 and 3, Supplementary Figure S1, Figure 1 data panels (ggplot2/ggseg)
```

## Input files

Place the files below in `data/ADNI_PVED/`, `data/HABS_PVED/` and `data/`, or point `config.py` to their locations (environment variables `PVED_ADNI_DIR`, `PVED_HABS_DIR`, `PVED_PHC_TAU_FILE`).

**ADNI** (`PVED_ADNI_DIR`)
- `ADNI_PVED_AandD_github_package_20260702/results/alz20_pacc_reanalysis/baseline_main_with_alz20_pacc.csv`: baseline PVeD analytic file with plasma biomarkers, demographics and the PACC composite, from the original submission's pipeline
- `pacc_long_pved.csv`: longitudinal PACC composite with time from the PVeD MRI
- `PVeD_all_merged.csv`: EstPVeD output
- `DTI_Image_Cohort_Analysis_Manifest_23Apr2026.csv`
- `UCBERKELEY_AMY_6MM_28Apr2026.csv`, `UCBERKELEY_TAU_6MM_28Apr2026.csv`
- `UCD_WMH_24May2026.csv`, `CDR_28Apr2026.csv`, `PTDEMOG_28Apr2026.csv`

**HABS-HD** (`PVED_HABS_DIR`)
- `PVeD_metrics.csv`: EstPVeD output
- `RP_HD_7_Clinical.csv`, `RP_HD_7_Biomarkers.xlsx`, `RP_HD_7_Genomics.xlsx`
- `HABS_HD_RP_8_FBB_PET_data_BraskieLab_20260501.xlsx`, `HABS_HD_RP_8_PI2620_PET_Data_BraskieLab_20260501.xlsx`
- `HABS_HD_RP_8_Hippocampal_Volume_HippoDeep_Data_BraskieLab_20260501.xlsx`, `HABS_HD_RP_8_CorticalThickness_ATN_MetaROI_Data_BraskieLab_20260501.xlsx`
- `HABS_HD_RP_8_WMH_Combination_Report_Analysis_Borzage_Lab_20260507.xlsx`, `HABS_HD_RP_8_Whole_Brain_Volume_Data_BraskieLab_20260501.xlsx`
- `HABS_HD_RP_8_DTI_ROI_Results.xlsx` (acquisition protocol only)
- `HABS_HD_ADNI_mirror_tables_figures_20260703/validation/HABS_pacc_like_cross_sectional_model_ready.csv` and `HABS_pacc_like_visit_level.csv`: PACC-like composite, from the original submission's pipeline

**ADSP-PHC** (`PVED_PHC_TAU_FILE`): `ADSP_PHC_PET_Tau_Detailed_*.csv`

File names carry the download dates used for the manuscript; adjust them in the scripts if your downloads differ.

## Running

```bash
pip install -r requirements.txt
python run_all.py                     # about 5 min; results/ holds one CSV per analysis
Rscript figures/R/make_figures.R      # needs R packages ggplot2, dplyr, readr, ggseg, scales
```

`python run_all.py --from lmm_profile.py` resumes from a given step.

## Outputs and where they appear in the manuscript

| File in `results/` | Manuscript |
|---|---|
| `T_table1.csv` | Table 1 |
| `T_models_M1_M4.csv` | Table 2, Table S4, Figure 3 |
| `T_lmm_continuous.csv`, `T_lmm_profiles.csv` | Table 3, Tables S14-S15, Figure 4 |
| `T_intervals.csv` | Table S3 |
| `T_missingness.csv` | Table S5 |
| `T_simple_slopes.csv` | Table S6, Supplementary Figure S2 |
| `T_time_windows.csv` | Table S7, Supplementary Figure S3 |
| `T_robustness.csv`, `T_extra_sensitivity.csv` | Tables S8, S9, S18 |
| `T_specificity.csv` | Table S10 |
| `T_incremental_R2.csv` | Table S11 |
| `T_amyloid_discrimination.csv`, `T_discordance.csv` | Tables S12-S13 |
| `T_two_stage_slopes.csv` | Table S16 |
| `T_followup.csv` | Table S17 |
| `T_tau_tracer_handling.csv` | Table S20 |
| `T_long_method_comparison.csv` | Table S21 |
| `T_long_meta.csv` | Table S22 |

Table S1 (cohort derivation counts) is compiled from the data-release manifests. Table S19 is unchanged from the original submission. Figure 1 is a schematic; only its data panels (sections 3 and 4) are generated here.

## Key analytic choices

- **Standardisation:** pTau217 (log10), PVeD, age and education are standardised within each cohort's analytic sample. The interaction term is the product of standardised pTau217 and PVeD.
- **PET timing:** PET measures are those closest in time to the PVeD MRI.
- **Tau PET tracers:**
  - ADNI tau PET analyses use flortaucipir only, because SUVRs are not comparable across tracers.
  - ADNI amyloid SUVR is standardised within tracer; Centiloid is the cross-tracer amyloid measure.
- **Covariates:**
  - Primary (M1): ADNI, age, sex and education; HABS-HD, age, sex and APOE ε4.
  - M2: harmonised covariates in both cohorts.
  - M3: M2 plus CDR sum of boxes.
- **Multiple comparisons:** Benjamini-Hochberg correction within cohort, model and endpoint family.
- **ADNI PACC:** mean of available z-scores for ADAS-Cog delayed word recall (reversed), MMSE and Logical Memory delayed recall, referenced to ADNI baseline cognitively normal participants.
- **HABS-HD PACC-like composite:** mean of internally standardised MMSE, SEVLT delayed recall, Digit Symbol Substitution and animal fluency.
- **Longitudinal models:** linear mixed-effects models (REML) with random intercepts and slopes for time and covariate-by-time terms. Outcomes are expressed in baseline-SD units per year.

## Contact

Wencai Ding, The Second Affiliated Hospital of Wannan Medical College (dingneurology@gmail.com).

## Revision-code entry point and earlier versions

The [`revision_analyses/`](revision_analyses/) entry indexes the resubmission scripts under `analysis/` and `figures/`. The earlier `R/`, `scripts/` and `src/` files remain for provenance; use `run_all.py` and the instructions above for the current resubmission workflow.
