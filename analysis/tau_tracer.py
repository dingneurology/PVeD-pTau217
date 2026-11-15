"""ADNI tau tracer handling: FTP-only primary vs all-tracer sensitivity analyses incl. ADSP-PHC normative z-scores."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import pandas as pd, numpy as np, sys, warnings
warnings.filterwarnings("ignore")
from statlib import *
from analyses import fit_int
R = CFG.ADNI
PHC = CFG.PHC_TAU
A = pd.read_pickle(CFG.WORK + "A.pkl")
d = pd.read_pickle(CFG.WORK + "adni.pkl")[["RID", "DTI_DATE"]]
t = pd.read_csv(R + "UCBERKELEY_TAU_6MM_28Apr2026.csv", low_memory=False); t["SCANDATE"] = pd.to_datetime(t.SCANDATE)
m = d.merge(t[["RID", "SCANDATE", "qc_flag"]], on="RID"); m["ad"] = (m.SCANDATE - m.DTI_DATE).dt.days.abs()
m = m.sort_values(["RID", "ad"]).drop_duplicates("RID"); m["id"] = m.RID.astype(str)
A = A.merge(m[["id", "qc_flag"]], on="id", how="left")
# ADSP-PHC (MRI-free pipeline; Level 2 SUVRs) -> tracer-specific normative z using PHC methods Table 7 (T- class mean, SD)
p = pd.read_csv(PHC, low_memory=False); p["SCANDATE"] = pd.to_datetime(p.PHC_SCANDATE)
mp = d.merge(p[p.PHC_QC_IMAGE == 1], on="RID"); mp["ad"] = (mp.SCANDATE - mp.DTI_DATE).dt.days.abs()
mp = mp.sort_values(["RID", "ad"]).drop_duplicates("RID"); mp["id"] = mp.RID.astype(str)
ref = {"FTP": {"B12": (1.15, .09), "B3": (1.17, .07), "B4": (1.17, .08), "B5": (1.10, .07), "B6": (1.02, .07)},
       "MK": {"B12": (1.01, .09), "B3": (1.03, .07), "B4": (1.02, .08), "B5": (0.98, .07), "B6": (0.91, .06)},
       "PI": {"B12": (1.13, .09), "B3": (1.10, .10), "B4": (1.04, .08), "B5": (0.97, .08), "B6": (0.95, .09)}}
for b in ref["FTP"]:
    mp["z" + b] = [(s - ref[tr][b][0]) / ref[tr][b][1] for s, tr in zip(mp[b + "_SUVR"], mp.PHC_TRACER)]
    mp["pos" + b] = (mp["z" + b] > 2.5).astype(float)
mp["zB34"] = mp[["zB3", "zB4"]].mean(axis=1)
def stage(r):
    s = 0
    for k, b in enumerate(["B12", "B3", "B4", "B5", "B6"], 1):
        if r["pos" + b] == 1 and s == k - 1:
            s = k
    return s
mp["hstage"] = mp.apply(stage, axis=1)
A = A.merge(mp[["id", "PHC_TRACER", "zB12", "zB3", "zB4", "zB5", "zB34", "hstage"]], on="id", how="left")
cov = ["za", "male", "ze"]
lab = {"tau_meta": "Meta-temporal SUVR", "tau_ent": "Entorhinal SUVR", "tau_it": "Inferior temporal SUVR", "tau_mt": "Middle temporal SUVR",
       "tau_fus": "Fusiform SUVR", "tau_pcun": "Precuneus SUVR", "tau_pcc": "Posterior cingulate SUVR",
       "zB12": "Braak I–II z", "zB3": "Braak III z", "zB4": "Braak IV z", "zB5": "Braak V z", "zB34": "Braak III–IV z", "hstage": "Hierarchical Braak stage (0–5)"}
taus = ["tau_meta", "tau_ent", "tau_it", "tau_mt", "tau_fus", "tau_pcun", "tau_pcc"]
rows = []
def add(name, df, ys, sub=None, cat=()):
    for y in ys:
        _, r = fit_int(df, y, cov, sub=sub, extra_cat=cat)
        rows.append(dict(analysis=name, outcome=lab[y], var=y, **r))
add("Primary: flortaucipir only (UC Berkeley SUVR)", A, taus)
add("Flortaucipir only, PET QC pass", A, taus, sub=A.qc_flag == 2)
B = A.copy()
for y in taus:
    B[y] = B[y + "_alltracer"]
add("All tracers pooled, unharmonized SUVR (original approach)", B, taus)
C = B.copy()
for y in taus:
    C[y] = C.groupby("tau_tracer")[y].transform(lambda s: (s - s.mean()) / s.std())
add("All tracers, within-tracer standardized SUVR + tracer covariate", C, taus, cat=("tau_tracer",))
add("All tracers, ADSP-PHC tracer-specific normative z (Braak ROIs)", A, ["zB12", "zB3", "zB4", "zB5", "zB34", "hstage"])
T = pd.DataFrame(rows)
am = [fit_int(A, y, cov)[1]["p"] for y in ["amy_cl", "amy_suvr"]]
T["fdr"] = np.nan
for k, g in T.groupby("analysis", sort=False):
    if g["var"].iloc[0].startswith("tau_"):
        T.loc[g.index, "fdr"] = bh(np.r_[am, g.p.values])[2:]
    else:
        T.loc[g.index, "fdr"] = bh(g.p.values)
T.to_csv(CFG.OUT + "T_tau_tracer_handling.csv", index=False)
cnt = dict(berkeley=A.tau_tracer.value_counts().to_dict(), phc=A.PHC_TRACER.value_counts().to_dict())
pd.Series({k: str(v) for k, v in cnt.items()}).to_csv(CFG.OUT + "T_tau_tracer_counts.csv")
print(cnt)
pd.set_option("display.width", 220)
print(T[["analysis", "outcome", "n", "beta", "lo", "hi", "p", "fdr"]].round(3).to_string())
