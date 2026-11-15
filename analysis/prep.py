"""Build harmonised analytic datasets for ADNI (discovery) and HABS-HD (validation)."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import re
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
from statlib import z  # noqa: E402

RA = CFG.ADNI
RH = CFG.HABS


def adni():
    d = pd.read_pickle(CFG.WORK + "adni.pkl").copy()
    d = d.copy()
    # QSDR registration goodness-of-fit (R) encoded in the DSI Studio file name
    pv = pd.read_csv(RA + "PVeD_all_merged.csv")
    pv["qsdr_R"] = pv.image_id.str.extract(r"\.R(\d+)\.fib")[0].astype(float)
    d = d.merge(pv[["RID", "qsdr_R"]].drop_duplicates("RID"), on="RID", how="left")
    # phase / protocol
    ptd = pd.read_csv(RA + "PTDEMOG_28Apr2026.csv", usecols=["RID", "PHASE", "SITEID", "VISDATE"]).sort_values(["RID", "VISDATE"]).drop_duplicates("RID")
    man = pd.read_csv(RA + "DTI_Image_Cohort_Analysis_Manifest_23Apr2026.csv", usecols=["image_id", "image_description"])
    d = d.merge(ptd[["RID", "PHASE"]], on="RID", how="left").merge(man, left_on="DTI_IMAGE_ID", right_on="image_id", how="left")
    out = pd.DataFrame({
        "id": d.RID.astype(str),
        "cohort": "ADNI",
        "age": d.AGE, "male": (d.SEX == 1).astype(float), "edu": d.EDUCATION,
        "apoe4": d.APOE4_carrier,
        "dx": d.DIAGNOSIS.map({"CN": 0, "MCI": 1, "Dementia": 2}),
        "cdrsb": d.cdr_CDRSB, "cdr_mo": d.cdr_dd / 30.4375,
        "race_white": (pd.to_numeric(d.RACE, errors="coerce") == 5).astype(float), "race_black": (pd.to_numeric(d.RACE, errors="coerce") == 4).astype(float),
        "hispanic": (pd.to_numeric(d.ETHNICITY, errors="coerce") == 1).astype(float),
        "ptau_raw": d.bio_plasma_ptau217, "pved": d.PVeD_total, "pved_L": d.PVeD_L, "pved_R": d.PVeD_R,
        "qa": d.QA_index, "qsdr_R": d.qsdr_R,
        "ptau_mo": d.bio_plasma_ptau217_days_from_dti / 30.4375,
        "amy_cl": d.amy_CENTILOIDS, "amy_suvr": d.amy_SUMMARY_SUVR,
        "amy_pos": d.amy_AMYLOID_STATUS, "amy_mo": d.amy_dd / 30.4375, "amy_tracer": d.amy_TRACER,
        "tau_meta": d.tau_META_TEMPORAL_SUVR, "tau_ent": d.tau_CTX_ENTORHINAL_SUVR,
        "tau_it": d.tau_CTX_INFERIORTEMPORAL_SUVR, "tau_mt": d.tau_CTX_MIDDLETEMPORAL_SUVR,
        "tau_fus": d.tau_CTX_FUSIFORM_SUVR, "tau_pcun": d.tau_CTX_PRECUNEUS_SUVR,
        "tau_pcc": d.tau_CTX_POSTERIORCINGULATE_SUVR, "tau_mo": d.tau_dd / 30.4375,
        "hippo": d.mri_TOTAL_HIPPO, "icv": d.mri_CEREBRUM_TCV, "wmh": np.log(d.mri_TOTAL_WMH + 1),
        "csf_frac": d.mri_TOTAL_CSF / d.mri_CEREBRUM_TCV if "mri_TOTAL_CSF" in d else np.nan,
        "brain_frac": d.mri_TOTAL_BRAIN / d.mri_CEREBRUM_TCV, "mri_mo": d.mri_dd / 30.4375,
        "pacc": d.PACC_COMPOSITE, "pacc_3comp": d.PACC_COMPOSITE.where(d.N_PACC_COMPONENTS >= 3), "pacc_ncomp": d.N_PACC_COMPONENTS, "mmse": d.MMSE, "pacc_mo": d.PACC_DAYS_FROM_DTI_SIGNED / 30.4375,
        "gfap_raw": d.bio_plasma_gfap, "nfl_raw": d.bio_plasma_nfl,
        "phase": d.PHASE, "protocol": d.image_description,
        "tau_tracer": d.tau_TRACER,
    })
    # Revision: tau PET SUVRs are not interchangeable across tracers ("DO NOT COMPARE SUVRs ACROSS TRACERS").
    # Primary ADNI tau analyses use flortaucipir (FTP) only; all-tracer copies kept for sensitivity analyses.
    taucols = ["tau_meta", "tau_ent", "tau_it", "tau_mt", "tau_fus", "tau_pcun", "tau_pcc", "tau_mo"]
    for c in taucols:
        out[c + "_alltracer"] = out[c]
        out.loc[out.tau_tracer != "FTP", c] = np.nan
    # Amyloid SUVR is tracer-specific: standardize within tracer before pooling (Centiloid is the cross-tracer measure)
    out["amy_suvr_raw"] = out.amy_suvr
    out["amy_suvr"] = out.groupby("amy_tracer").amy_suvr_raw.transform(lambda s: (s - s.mean()) / s.std())
    return out


def habs():
    d = pd.read_pickle(CFG.WORK + "habs_raw.pkl")
    d = d[d.ptau_plus.notna()].copy()
    cl = pd.read_csv(RH + "RP_HD_7_Clinical.csv", low_memory=False, usecols=lambda c: c in ["Med_ID", "Visit_ID", "ID_Race_White", "ID_Race_Black", "ID_Hispanic", "ID_Education"])
    cl["Med_ID"] = cl.Med_ID.astype(str).str.extract(r"^(\d+)", expand=False)
    cl = cl.sort_values(["Med_ID", "Visit_ID"]).drop_duplicates("Med_ID")
    d = d.drop(columns=[c for c in ["ID_Hispanic", "ID_Education"] if c in d]).merge(cl.drop(columns="Visit_ID"), on="Med_ID", how="left")
    for c in ["ID_Race_White", "ID_Race_Black", "ID_Hispanic", "ID_Education", "APOE4_Positivity", "CDR_Sum_pv", "MMSE_Total", "CDX_Cog_pv", "ID_Gender"]:
        d[c] = pd.to_numeric(d[c], errors="coerce").where(lambda x: x >= 0)
    # whole-brain volume (BrainSegNotVent / ICV) nearest to the PVeD MRI
    xl = pd.ExcelFile(RH + "HABS_HD_RP_8_Whole_Brain_Volume_Data_BraskieLab_20260501.xlsx")
    wb = pd.concat([pd.read_excel(xl, sheet_name=s) for s in xl.sheet_names if "Included" in s])
    wb["Med_ID"] = wb.Med_ID.astype(str).str.extract(r"^(\d+)", expand=False)
    wb["wb_date"] = pd.to_datetime(dict(year=pd.to_numeric(wb.MRI_year, errors="coerce"), month=pd.to_numeric(wb.MRI_month, errors="coerce"), day=15), errors="coerce")
    wb = wb.merge(d[["Med_ID", "pved_mri_date"]], on="Med_ID")
    wb["amo"] = (wb.wb_date - wb.pved_mri_date).dt.days.abs()
    wb = wb.sort_values(["Med_ID", "amo"]).drop_duplicates("Med_ID")
    wb["brain_frac"] = pd.to_numeric(wb.BrainSegNotVent, errors="coerce") / pd.to_numeric(wb.ICV_HD, errors="coerce")
    d = d.merge(wb[["Med_ID", "brain_frac"]], on="Med_ID", how="left")
    # DTI protocol nearest to PVeD MRI
    dti = pd.read_pickle(CFG.WORK + "habs_dti.pkl")
    dti["Med_ID"] = dti.Med_ID.astype(str).str.extract(r"^(\d+)", expand=False)
    dti["ddate"] = pd.to_datetime(dict(year=dti.Year, month=dti.Month, day=15), errors="coerce")
    dti = dti.merge(d[["Med_ID", "pved_mri_date"]], on="Med_ID")
    dti["amo"] = (dti.ddate - dti.pved_mri_date).dt.days.abs()
    dti = dti.sort_values(["Med_ID", "amo"]).drop_duplicates("Med_ID")
    d = d.merge(dti[["Med_ID", "Protocol"]], on="Med_ID", how="left")
    # cross-sectional PACC-like composite (closest visit to PVeD MRI, released mirror file)
    mr = pd.read_csv(RH + "HABS_HD_ADNI_mirror_tables_figures_20260703/validation/HABS_pacc_like_cross_sectional_model_ready.csv", low_memory=False)
    mr["Med_ID"] = mr.Med_ID.astype(str).str.extract(r"^(\d+)", expand=False)
    mr = mr.drop_duplicates("Med_ID")[["Med_ID", "HABS_PACC_like_primary", "HABS_PACC_like_strict", "HABS_PACC_like_trails_sensitivity", "HABS_PACC_like_neuropsych_sensitivity", "pacc_months_from_pved"]]
    d = d.merge(mr, on="Med_ID", how="left")
    hip = d.H_L_hippocampus + d.H_R_hippocampus
    out = pd.DataFrame({
        "id": d.Med_ID.astype(str), "cohort": "HABS-HD",
        "age": d.Age, "male": (d.ID_Gender == 0).astype(float).where(d.ID_Gender.notna()), "edu": d.ID_Education,
        "apoe4": d.APOE4_Positivity, "dx": d.CDX_Cog_pv.where(d.CDX_Cog_pv <= 2),
        "cdrsb": d.CDR_Sum_pv, "cdr_mo": d.ptau_mo,
        "race_white": d.ID_Race_White, "race_black": d.ID_Race_Black,
        "hispanic": (d.ID_Hispanic >= 2).astype(float).where(d.ID_Hispanic.notna()),
        "ptau_raw": d.ptau_plus, "pved": d.PVeD_total, "pved_L": d.PVeD_L, "pved_R": d.PVeD_R, "qa": d.QA_index,
        "ptau_mo": d.ptau_mo,
        "amy_cl": d.Anear_Centiloid_Units_GlobalAmyloid, "amy_suvr": d.Anear_Global_FBB_SUVR,
        "amy_pos": d.Anear_AB_pos_FBB, "amy_mo": d.Anear_mo,
        "tau_mtl": d.Tnear_Medial_Temporal_PI2620_SUVR, "tau_pcc": d.Tnear_Posterior_Cingulate_PI2620_SUVR,
        "tau_lp": d.Tnear_Lateral_Parietal_PI2620_SUVR, "tau_mstl": d.Tnear_Tennant_MsTL_taucutpoint_weighted_SUVR,
        "tau_mo": d.Tnear_mo,
        "hippo": hip, "icv": d.H_ICV_HD, "ct_meta": d.CT_CorticalThickness_MetaROI, "wmh": d.W_wmh_log,
        "brain_frac": d.brain_frac,
        "pacc": d.HABS_PACC_like_primary, "pacc_strict": d.HABS_PACC_like_strict,
        "pacc_trails": d.HABS_PACC_like_trails_sensitivity, "pacc_np": d.HABS_PACC_like_neuropsych_sensitivity,
        "pacc_mo": d.pacc_months_from_pved, "mmse": pd.to_numeric(d.MMSE_Total, errors="coerce").where(lambda x: x >= 0),
        "gfap_raw": d.gfap_plus, "nfl_raw": d.nfl_r3,
        "gfap_same_visit": (d.gfap_plus_visit == d.ptau_plus_visit).astype(float),
        "nfl_same_visit": (d.nfl_r3_visit == d.ptau_plus_visit).astype(float),
        "protocol": d.Protocol, "pved_mri_date": d.pved_mri_date,
    })
    return out


def finish(df):
    df = df.copy()
    df["ptau_log"] = np.log10(df.ptau_raw.where(df.ptau_raw > 0))
    df["zp"] = z(df.ptau_log)
    df["zv"] = z(df.pved)
    df["za"] = z(df.age)
    df["ze"] = z(df.edu)
    df["int"] = df.zp * df.zv
    df["zcdr"] = z(df.cdrsb)
    df["zicv"] = z(df.icv)
    for m in ["gfap", "nfl"]:
        v = df[m + "_raw"]
        df["z" + m[0]] = z(np.log10(v.where(v > 0)))
    df["zg_int"] = df.zg * df.zv
    df["zn_int"] = df.zn * df.zv
    df["hi_ptau"] = (df.ptau_log > df.ptau_log.median()).astype(int)
    df["hi_pved"] = (df.pved > df.pved.median()).astype(int)
    lab = {(0, 1): "Low pTau217 / High PVeD", (0, 0): "Low pTau217 / Low PVeD", (1, 1): "High pTau217 / High PVeD", (1, 0): "High pTau217 / Low PVeD"}
    df["profile"] = [lab[(a, b)] for a, b in zip(df.hi_ptau, df.hi_pved)]
    df["impaired"] = (df.dx >= 1).astype(float).where(df.dx.notna())
    return df


if __name__ == "__main__":
    A = finish(adni())
    H = finish(habs())
    A.to_pickle(CFG.WORK + "A.pkl")
    H.to_pickle(CFG.WORK + "H.pkl")
    print(A.shape, H.shape)
    print(A.profile.value_counts().to_dict())
    print(H.profile.value_counts().to_dict())
    print("ADNI medians", A.ptau_raw.median(), A.pved.median(), "HABS", H.ptau_raw.median(), H.pved.median())
    print(A[["cdrsb", "apoe4", "race_white", "hispanic", "csf_frac", "qsdr_R"]].describe().T)
    print(H[["cdrsb", "apoe4", "edu", "race_white", "race_black", "hispanic", "brain_frac", "gfap_raw", "nfl_raw"]].describe().T)
