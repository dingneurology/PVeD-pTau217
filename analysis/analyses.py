"""Revision analyses for BMC Medicine manuscript (PVeD x pTau217).
All outputs are written as CSV files to the results/ folder (see config.py).
"""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
from statlib import (auc, bh, delong, design, lincom, lmm_fit, logit, nri_idi, ols, z)  # noqa: E402

OUT = CFG.OUT
A = pd.read_pickle(CFG.WORK + "A.pkl")
H = pd.read_pickle(CFG.WORK + "H.pkl")

OUTC = {
    "ADNI": [
        ("Amyloid PET", "amy_cl", "Amyloid PET Centiloid"),
        ("Amyloid PET", "amy_suvr", "Amyloid PET SUVR"),
        ("Tau PET", "tau_meta", "Tau PET meta-temporal SUVR"),
        ("Tau PET", "tau_ent", "Tau PET entorhinal SUVR"),
        ("Tau PET", "tau_it", "Tau PET inferior temporal SUVR"),
        ("Tau PET", "tau_mt", "Tau PET middle temporal SUVR"),
        ("Tau PET", "tau_fus", "Tau PET fusiform SUVR"),
        ("Tau PET", "tau_pcun", "Tau PET precuneus SUVR"),
        ("Tau PET", "tau_pcc", "Tau PET posterior cingulate SUVR"),
        ("MRI", "hippo", "Hippocampal volume"),
        ("MRI", "wmh", "WMH volume (log)"),
        ("Cognition", "pacc", "Baseline PACC"),
        ("Cognition", "mmse", "Baseline MMSE"),
    ],
    "HABS-HD": [
        ("Amyloid PET", "amy_cl", "FBB amyloid Centiloid"),
        ("Amyloid PET", "amy_suvr", "FBB global amyloid SUVR"),
        ("Tau PET", "tau_mtl", "PI-2620 medial temporal SUVR"),
        ("Tau PET", "tau_pcc", "PI-2620 posterior cingulate SUVR"),
        ("Tau PET", "tau_lp", "PI-2620 lateral parietal SUVR"),
        ("Tau PET", "tau_mstl", "PI-2620 MsTL tau-stage weighted SUVR"),
        ("MRI", "hippo", "Hippocampal volume"),
        ("MRI", "ct_meta", "AD-signature cortical thickness"),
        ("MRI", "wmh", "WMH volume (log)"),
        ("Cognition", "pacc", "Baseline PACC"),
        ("Cognition", "mmse", "Baseline MMSE"),
    ],
}
DATA = {"ADNI": A, "HABS-HD": H}
PRIMARY_COV = {"ADNI": ["za", "male", "ze"], "HABS-HD": ["za", "male", "apoe4"]}
HARM_COV = ["za", "male", "ze", "apoe4"]


def fmt(b, lo, hi, d=3):
    return f"{b:.{d}f} ({lo:.{d}f} to {hi:.{d}f})"


def fit_int(df, y, covs, term="int", main=("zp", "zv"), extra_cat=(), sub=None):
    dd = df if sub is None else df[sub]
    dd = dd.copy()
    dd["_y"] = z(dd[y])
    cv = list(covs)
    if y == "hippo":
        cv = cv + ["zicv"]
    X = design(dd, list(main) + [term] + cv, cat=extra_cat)
    f = ols(dd["_y"], X)
    r = f["coef"].loc[term]
    return f, dict(n=f["n"], beta=r.beta, lo=r.lo, hi=r.hi, p=r.p, se=r.se)


# ---------------------------------------------------------------- 1. primary + harmonised + CDR models
def models_table():
    rows = []
    for coh, lst in OUTC.items():
        df = DATA[coh]
        sets = {
            "M1 primary (cohort-specific)": PRIMARY_COV[coh],
            "M2 harmonized (age, sex, education, APOE e4)": HARM_COV,
            "M3 M2 + CDR-SB": HARM_COV + ["zcdr"],
            "M4 M2 + clinical diagnosis": HARM_COV,
        }
        for mname, covs in sets.items():
            for dom, y, lab in lst:
                if mname.startswith("M3") and y in ("pacc", "mmse"):
                    pass
                cat = ("dx",) if mname.startswith("M4") else ()
                _, r = fit_int(df, y, covs, extra_cat=cat)
                rows.append(dict(cohort=coh, model=mname, domain=dom, outcome=lab, var=y, **r))
    t = pd.DataFrame(rows)
    t["family"] = np.where(t.domain.isin(["Amyloid PET", "Tau PET"]), "Primary (PET)", "Secondary (MRI/cognition)")
    t["fdr"] = np.nan
    t["fdr_all"] = np.nan
    for (c, m, fam), g in t.groupby(["cohort", "model", "family"]):
        t.loc[g.index, "fdr"] = bh(g.p.values)
    for (c, m), g in t.groupby(["cohort", "model"]):
        t.loc[g.index, "fdr_all"] = bh(g.p.values)
    t["est"] = [fmt(*x) for x in t[["beta", "lo", "hi"]].values]
    t.to_csv(OUT + "T_models_M1_M4.csv", index=False)
    return t


# ---------------------------------------------------------------- 2. simple slopes
def simple_slopes():
    rows = []
    for coh, lst in OUTC.items():
        df = DATA[coh]
        for dom, y, lab in lst:
            if dom not in ("Amyloid PET", "Tau PET"):
                continue
            f, r = fit_int(df, y, PRIMARY_COV[coh])
            out = dict(cohort=coh, outcome=lab, n=f["n"], int_beta=r["beta"])
            for name, v in [("low", -1), ("mean", 0), ("high", 1)]:
                e, lo, hi, p = lincom(f, {"zp": 1, "int": v})
                out.update({f"slope_{name}": e, f"slope_{name}_lo": lo, f"slope_{name}_hi": hi, f"slope_{name}_p": p})
            out["ratio_low_high"] = out["slope_low"] / out["slope_high"]
            rows.append(out)
    t = pd.DataFrame(rows)
    t.to_csv(OUT + "T_simple_slopes.csv", index=False)
    return t


# ---------------------------------------------------------------- 3. temporal alignment
def temporal():
    rows = []
    for coh, df in DATA.items():
        for var, lab in [("ptau_mo", "Plasma pTau217 to PVeD MRI"), ("amy_mo", "Amyloid PET to PVeD MRI"), ("tau_mo", "Tau PET to PVeD MRI"), ("pacc_mo", "Baseline cognition to PVeD MRI")]:
            x = df[var].dropna()
            ax = x.abs()
            rows.append(dict(cohort=coh, interval=lab, n=len(x), median_abs_months=ax.median(), q1=ax.quantile(.25), q3=ax.quantile(.75), max_abs=ax.max(),
                             signed_median=x.median(), pct_within_6=(ax <= 6).mean() * 100, pct_within_12=(ax <= 12).mean() * 100, pct_within_24=(ax <= 24).mean() * 100))
    pd.DataFrame(rows).to_csv(OUT + "T_intervals.csv", index=False)
    rows = []
    wins = {"ADNI": [None, 12, 6], "HABS-HD": [None, 24, 12, 6]}
    for coh, lst in OUTC.items():
        df = DATA[coh]
        for w in wins[coh]:
            for dom, y, lab in lst:
                if dom not in ("Amyloid PET", "Tau PET"):
                    continue
                mo = "amy_mo" if dom == "Amyloid PET" else "tau_mo"
                sub = None if w is None else ((df[mo].abs() <= w) & (df.ptau_mo.abs() <= w))
                _, r = fit_int(df, y, PRIMARY_COV[coh], sub=sub)
                rows.append(dict(cohort=coh, window="All available" if w is None else f"±{w} months", outcome=lab, **r))
            # interval as covariate
        for dom, y, lab in lst:
            if dom not in ("Amyloid PET", "Tau PET"):
                continue
            mo = "amy_mo" if dom == "Amyloid PET" else "tau_mo"
            dd = df.copy(); dd["abs_gap"] = dd[mo].abs(); dd["abs_pgap"] = dd.ptau_mo.abs()
            _, r = fit_int(dd, y, PRIMARY_COV[coh] + ["abs_gap", "abs_pgap"])
            rows.append(dict(cohort=coh, window="All available + interval covariates", outcome=lab, **r))
    t = pd.DataFrame(rows)
    t["fdr"] = np.nan
    for k, g in t.groupby(["cohort", "window"]):
        t.loc[g.index, "fdr"] = bh(g.p.values)
    t.to_csv(OUT + "T_time_windows.csv", index=False)
    return t


# ---------------------------------------------------------------- 4. longitudinal
def load_long():
    RA = CFG.ADNI
    L = pd.read_csv(RA + "pacc_long_pved.csv")
    L["id"] = L.RID.astype(str)
    # ADNI PACC: three-test composite analysed in the original submission (mean of available z-scores for ADAS-Cog delayed
    # word recall [reversed], MMSE and Logical Memory delayed recall, referenced to ADNI baseline CN); baseline-SD units
    L = L[L.id.isin(A.id)][["id", "time_from_dti_years", "PACC_COMPOSITE"]].rename(columns={"time_from_dti_years": "t", "PACC_COMPOSITE": "y"}).dropna(subset=["t", "y"])
    L["y"] = L.y / A.pacc.std()
    LA = L.merge(A, on="id")
    RH = CFG.HABS
    v = pd.read_csv(RH + "HABS_HD_ADNI_mirror_tables_figures_20260703/validation/HABS_pacc_like_visit_level.csv", low_memory=False)
    v["id"] = v.Med_ID.astype(str).str.extract(r"^(\d+)", expand=False)
    v = v.merge(H[["id", "pved_mri_date"]], on="id")
    v["t"] = (pd.to_datetime(v.Visit_Date_parsed) - v.pved_mri_date).dt.days / 365.25
    v = v[["id", "t", "HABS_PACC_like_primary", "Visit_Date_parsed"]].rename(columns={"HABS_PACC_like_primary": "y"}).dropna(subset=["y", "t"])
    v["y"] = v.y / H.pacc.std()
    LH = v.merge(H.drop(columns=["pved_mri_date"]), on="id")
    return LA, LH


def lmm_three_way(L, covs, label, min_visits=1):
    L = L.copy()
    cnt = L.groupby("id").size()
    L = L[L.id.isin(cnt[cnt >= min_visits].index)]
    terms = ["zp", "zv", "int"] + covs
    L = L.dropna(subset=["y", "t"] + terms)
    # standardise outcome by baseline-visit SD so slopes are in baseline-SD units/year? keep native PACC units
    Xd = {"const": np.ones(len(L)), "t": L.t.values}
    for c in terms:
        Xd[c] = L[c].values.astype(float)
        Xd[c + ":t"] = L[c].values.astype(float) * L.t.values
    X = pd.DataFrame(Xd)
    Zm = np.column_stack([np.ones(len(L)), L.t.values])
    f = lmm_fit(L.y.values, X.values, Zm, L.id.values)
    names = list(X.columns)
    cov = pd.DataFrame(f["cov"], index=names, columns=names)
    b = pd.Series(f["beta"], index=names)
    dfree = len(L) - len(names)

    def lc(w):
        v = np.array([w.get(k, 0) for k in names])
        e = float(v @ b.values); se = float(np.sqrt(v @ cov.values @ v))
        return e, e - 1.96 * se, e + 1.96 * se, 2 * stats.norm.sf(abs(e / se))

    res = dict(model=label, n_participants=L.id.nunique(), n_obs=len(L), median_fu=L.groupby("id").t.agg(lambda x: x.max() - x.min()).median(), converged=f["converged"])
    e, lo, hi, p = lc({"int:t": 1})
    res.update(three_way=e, three_way_lo=lo, three_way_hi=hi, three_way_p=p)
    e, lo, hi, p = lc({"zp:t": 1})
    res.update(ptau_t=e, ptau_t_lo=lo, ptau_t_hi=hi, ptau_t_p=p)
    # annual change at +1 SD pTau217 for PVeD -1 vs +1 SD: difference = (zv:t)*(-2) + (int:t)*(1*-2)
    e, lo, hi, p = lc({"zv:t": -2, "int:t": -2})
    res.update(diff_hi_ptau_lowVShigh_pved=e, diff_lo=lo, diff_hi=hi, diff_p=p)
    # baseline (intercept) interaction
    e, lo, hi, p = lc({"int": 1})
    res.update(int_baseline=e, int_baseline_lo=lo, int_baseline_hi=hi, int_baseline_p=p)
    # predicted slopes in 4 corners
    for pl, pv_ in [(1, -1), (1, 1), (-1, -1), (-1, 1)]:
        e, lo, hi, p = lc({"t": 1, "zp:t": pl, "zv:t": pv_, "int:t": pl * pv_})
        res[f"slope_ptau{pl:+d}_pved{pv_:+d}"] = e
        res[f"slope_ptau{pl:+d}_pved{pv_:+d}_ci"] = f"{lo:.3f} to {hi:.3f}"
    return res


def longitudinal():
    LA, LH = load_long()
    rows = []
    for coh, L, base in [("ADNI", LA, ["za", "male", "ze"]), ("HABS-HD", LH, ["za", "male", "ze"])]:
        specs = [
            ("LMM-1: age, sex, education", base, L, 1),
            ("LMM-2: + APOE e4", base + ["apoe4"], L, 1),
            ("LMM-3: + APOE e4 + baseline CDR-SB", base + ["apoe4", "zcdr"], L, 1),
            ("LMM-4: participants with >=2 assessments", base, L, 2),
        ]
        if coh == "HABS-HD":
            post = L[L.t >= -0.5].copy()
            specs.append(("LMM-5: assessments from PVeD MRI onward only (t >= -6 months)", base, post, 1))
        else:
            specs.append(("LMM-5: assessments from PVeD MRI onward only (all ADNI visits already t >= 0)", base, L[L.t >= 0], 1))
        for lab, cv, dat, mv in specs:
            r = lmm_three_way(dat, cv, lab, mv)
            r["cohort"] = coh
            rows.append(r)
            print(coh, lab, round(r["three_way"], 4), round(r["three_way_p"], 4), r["n_participants"], flush=True)
    t = pd.DataFrame(rows)
    t["fdr_three_way"] = np.nan
    for c, g in t.groupby("cohort"):
        t.loc[g.index, "fdr_three_way"] = bh(g.three_way_p.values)
    t.to_csv(OUT + "T_lmm_continuous.csv", index=False)
    # follow-up structure
    fu = []
    for coh, L, base_df in [("ADNI", LA, A), ("HABS-HD", LH, H)]:
        Lp = L if coh == "ADNI" else L[L.t >= -0.5]
        g = Lp.groupby("id").agg(n=("t", "size"), fu=("t", lambda x: x.max() - x.min()))
        g2 = g[g.n >= 2]
        fu.append(dict(cohort=coh, baseline_n=len(base_df), with_any_assessment=len(g), with_2plus=len(g2), pct_2plus=len(g2) / len(base_df) * 100,
                       visits_mean=g2.n.mean(), visits_median=g2.n.median(), visits_max=g2.n.max(), fu_mean=g2.fu.mean(), fu_sd=g2.fu.std(), fu_median=g2.fu.median(),
                       fu_q1=g2.fu.quantile(.25), fu_q3=g2.fu.quantile(.75), single_assessment=(g.n == 1).sum()))
        # by profile
        gp = g2.join(base_df.set_index("id").profile)
        for pf, gg in gp.groupby("profile"):
            fu.append(dict(cohort=coh, profile=pf, with_2plus=len(gg), visits_mean=gg.n.mean(), fu_mean=gg.fu.mean(), fu_median=gg.fu.median(), obs=gg.n.sum()))
    pd.DataFrame(fu).to_csv(OUT + "T_followup.csv", index=False)
    # two-stage profile slopes with harmonised covariates, and continuous slope regression
    rows = []
    for coh, L in [("ADNI", LA), ("HABS-HD", LH if True else None)]:
        Lp = L if coh == "ADNI" else L
        sl = []
        for i, g in Lp.groupby("id"):
            if len(g) >= 2 and g.t.max() - g.t.min() > 0:
                sl.append((i, stats.linregress(g.t, g.y).slope))
        S = pd.DataFrame(sl, columns=["id", "slope"]).merge(DATA[coh], on="id")
        for lab, cv in [("age, sex, education", ["za", "male", "ze"]), ("+ APOE e4", ["za", "male", "ze", "apoe4"]), ("+ APOE e4 + CDR-SB", ["za", "male", "ze", "apoe4", "zcdr"])]:
            S2 = S.copy()
            for pf in ["Low pTau217 / Low PVeD", "High pTau217 / High PVeD", "High pTau217 / Low PVeD"]:
                S2["p_" + pf] = (S2.profile == pf).astype(float)
            X = design(S2, ["p_Low pTau217 / Low PVeD", "p_High pTau217 / High PVeD", "p_High pTau217 / Low PVeD"] + cv)
            f = ols(S2.slope, X)
            e, lo, hi, p = lincom(f, {"p_High pTau217 / Low PVeD": 1, "p_High pTau217 / High PVeD": -1})
            # continuous
            X2 = design(S2, ["zp", "zv", "int"] + cv)
            f2 = ols(S2.slope, X2)
            r2 = f2["coef"].loc["int"]
            rows.append(dict(cohort=coh, covariates=lab, n=f["n"], contrast=e, lo=lo, hi=hi, p=p, cont_int=r2.beta, cont_lo=r2.lo, cont_hi=r2.hi, cont_p=r2.p))
    pd.DataFrame(rows).to_csv(OUT + "T_two_stage_slopes.csv", index=False)


# ---------------------------------------------------------------- 5. amyloid PET positivity / discordance
def cv_auc(df, ycol, Xcols, k=10, reps=50, seed=7):
    d = df.dropna(subset=[ycol] + Xcols)
    rng = np.random.default_rng(seed)
    aucs = []
    for r in range(reps):
        folds = rng.permutation(len(d)) % k
        pred = np.zeros(len(d))
        for f in range(k):
            tr, te = d[folds != f], d[folds == f]
            m = logit(tr[ycol], design(tr, Xcols))
            Xte = design(te, Xcols)[m["coef"].index].to_numpy(float)
            pred[folds == f] = 1 / (1 + np.exp(-(Xte @ m["coef"].beta.values)))
        aucs.append(auc(d[ycol].values, pred))
    return np.mean(aucs), np.percentile(aucs, 2.5), np.percentile(aucs, 97.5)


def discrimination():
    rows, disc = [], []
    for coh, df in DATA.items():
        d = df.copy()
        d["apos"] = d.amy_pos.astype(float)
        cov = ["za", "male", "ze", "apoe4"]
        need = ["apos", "zp", "zv"] + cov
        d = d.dropna(subset=need)
        specs = {"Covariates": cov, "+ pTau217": cov + ["zp"], "+ pTau217 + PVeD": cov + ["zp", "zv"], "+ pTau217 x PVeD": cov + ["zp", "zv", "int"],
                 "pTau217 only": ["zp"], "pTau217 x PVeD only": ["zp", "zv", "int"]}
        fits = {k: logit(d.apos, design(d, v)) for k, v in specs.items()}
        base = fits["+ pTau217"]
        for k, f in fits.items():
            a = auc(f["y"], f["prob"])
            cva = cv_auc(d, "apos", specs[k])
            rows.append(dict(cohort=coh, model=k, n=f["n"], n_pos=int(f["y"].sum()), auc=a, cv_auc=cva[0], cv_auc_lo=cva[1], cv_auc_hi=cva[2], loglik=f["ll"], aic=-2 * f["ll"] + 2 * f["k"]))
        dl = delong(base["y"], base["prob"], fits["+ pTau217 x PVeD"]["prob"])
        lr = 2 * (fits["+ pTau217 x PVeD"]["ll"] - base["ll"])
        nri, nri_e, nri_n, idi = nri_idi(base["y"], base["prob"], fits["+ pTau217 x PVeD"]["prob"])
        intc = fits["+ pTau217 x PVeD"]["coef"].loc["int"]
        rows.append(dict(cohort=coh, model="Comparison: + pTau217 x PVeD vs + pTau217", n=base["n"], auc=dl["auc1"], auc2=dl["auc2"], delta_auc=dl["diff"], delta_lo=dl["diff_lo"], delta_hi=dl["diff_hi"], delong_p=dl["p"],
                         lr_chi2=lr, lr_df=2, lr_p=stats.chi2.sf(lr, 2), nri=nri, nri_events=nri_e, nri_nonevents=nri_n, idi=idi,
                         logit_int_or=np.exp(intc.beta), logit_int_lo=np.exp(intc.beta - 1.96 * intc.se), logit_int_hi=np.exp(intc.beta + 1.96 * intc.se), logit_int_p=intc.p))
        # two-cutoff pTau217 approach (95% sensitivity lower, 95% specificity upper) on raw log pTau217
        x = d.ptau_log.values; y = d.apos.values
        lower = np.quantile(x[y == 1], 0.05)
        upper = np.quantile(x[y == 0], 0.95)
        d["zone"] = np.where(d.ptau_log < lower, "low", np.where(d.ptau_log > upper, "high", "intermediate"))
        for zone, g in d.groupby("zone"):
            rec = dict(cohort=coh, zone=zone, n=len(g), apos_n=int(g.apos.sum()), apos_pct=g.apos.mean() * 100, lower_cut=10 ** lower, upper_cut=10 ** upper)
            if zone == "intermediate" and g.apos.nunique() == 2 and len(g) > 40:
                f1 = logit(g.apos, design(g, cov + ["zp"]))
                f2 = logit(g.apos, design(g, cov + ["zp", "zv", "int"]))
                dl2 = delong(f1["y"], f1["prob"], f2["prob"])
                rec.update(auc_ptau=dl2["auc1"], auc_ptau_pved=dl2["auc2"], delta=dl2["diff"], delta_lo=dl2["diff_lo"], delta_hi=dl2["diff_hi"], p=dl2["p"])
                c1 = cv_auc(g, "apos", cov + ["zp"]); c2 = cv_auc(g, "apos", cov + ["zp", "zv", "int"])
                rec.update(cv_auc_ptau=c1[0], cv_auc_ptau_pved=c2[0])
                # PVeD main effect within intermediate zone
                f3 = logit(g.apos, design(g, cov + ["zp", "zv"]))
                pv = f3["coef"].loc["zv"]
                rec.update(pved_or_per_sd=np.exp(pv.beta), pved_or_lo=np.exp(pv.beta - 1.96 * pv.se), pved_or_hi=np.exp(pv.beta + 1.96 * pv.se), pved_or_p=pv.p)
            disc.append(rec)
        # plasma-high / PET-negative discordance by PVeD among those above median pTau217
        hp = d[d.hi_ptau == 1].copy()
        hp["discord"] = 1 - hp.apos
        f = logit(hp.discord, design(hp, cov + ["zp", "zv"]))
        pv = f["coef"].loc["zv"]
        disc.append(dict(cohort=coh, zone="High pTau217 (>median): odds of amyloid-PET-negative discordance per +1 SD PVeD", n=len(hp), apos_n=int(hp.apos.sum()),
                         apos_pct=hp.apos.mean() * 100, pved_or_per_sd=np.exp(pv.beta), pved_or_lo=np.exp(pv.beta - 1.96 * pv.se), pved_or_hi=np.exp(pv.beta + 1.96 * pv.se), pved_or_p=pv.p,
                         disc_pct_low_pved=(1 - hp[hp.hi_pved == 0].apos.mean()) * 100, disc_pct_high_pved=(1 - hp[hp.hi_pved == 1].apos.mean()) * 100))
        lp = d[d.hi_ptau == 0].copy()
        f = logit(lp.apos, design(lp, cov + ["zp", "zv"]))
        pv = f["coef"].loc["zv"]
        disc.append(dict(cohort=coh, zone="Low pTau217 (<=median): odds of amyloid-PET-positive discordance per +1 SD PVeD", n=len(lp), apos_n=int(lp.apos.sum()),
                         apos_pct=lp.apos.mean() * 100, pved_or_per_sd=np.exp(pv.beta), pved_or_lo=np.exp(pv.beta - 1.96 * pv.se), pved_or_hi=np.exp(pv.beta + 1.96 * pv.se), pved_or_p=pv.p,
                         disc_pct_low_pved=lp[lp.hi_pved == 0].apos.mean() * 100, disc_pct_high_pved=lp[lp.hi_pved == 1].apos.mean() * 100))
    pd.DataFrame(rows).to_csv(OUT + "T_amyloid_discrimination.csv", index=False)
    pd.DataFrame(disc).to_csv(OUT + "T_discordance.csv", index=False)


# ---------------------------------------------------------------- 6. biomarker specificity in identical participants
def specificity():
    rows = []
    for coh, df in DATA.items():
        lst = [x for x in OUTC[coh] if x[0] in ("Amyloid PET", "Tau PET")]
        d = df.copy()
        if coh == "HABS-HD":
            d = d[(d.gfap_same_visit == 1) | (d.nfl_same_visit == 1) | True]
        for dom, y, lab in lst:
            common = d.dropna(subset=[y, "zp", "zg", "zn", "zv"] + PRIMARY_COV[coh])
            common = common.copy()
            # re-standardise within the common sample for comparability
            for v in ["zp", "zg", "zn", "zv"]:
                common[v] = z(common[v])
            common["int"] = common.zp * common.zv; common["zg_int"] = common.zg * common.zv; common["zn_int"] = common.zn * common.zv
            for mlab, main, term in [("pTau217 x PVeD", ("zp", "zv"), "int"), ("GFAP x PVeD", ("zg", "zv"), "zg_int"), ("NfL x PVeD", ("zn", "zv"), "zn_int")]:
                _, r = fit_int(common, y, PRIMARY_COV[coh], term=term, main=main)
                rows.append(dict(cohort=coh, outcome=lab, model=mlab, term=mlab, **r))
            # joint model: all three biomarkers and all three interactions
            cc = common.copy(); cc["_y"] = z(cc[y])
            X = design(cc, ["zp", "zg", "zn", "zv", "int", "zg_int", "zn_int"] + PRIMARY_COV[coh])
            f = ols(cc["_y"], X)
            for term, tl in [("int", "pTau217 x PVeD"), ("zg_int", "GFAP x PVeD"), ("zn_int", "NfL x PVeD")]:
                r = f["coef"].loc[term]
                rows.append(dict(cohort=coh, outcome=lab, model="Joint model (all three biomarkers and interactions)", term=tl, n=f["n"], beta=r.beta, lo=r.lo, hi=r.hi, p=r.p, se=r.se))
            # pTau217 x PVeD + GFAP, NfL main effects only
            X = design(cc, ["zp", "zg", "zn", "zv", "int"] + PRIMARY_COV[coh])
            f = ols(cc["_y"], X); r = f["coef"].loc["int"]
            rows.append(dict(cohort=coh, outcome=lab, model="pTau217 x PVeD + GFAP and NfL main effects", term="pTau217 x PVeD", n=f["n"], beta=r.beta, lo=r.lo, hi=r.hi, p=r.p, se=r.se))
            # difference pTau217xPVeD - GFAPxPVeD in joint model
            Xj = design(cc, ["zp", "zg", "zn", "zv", "int", "zg_int", "zn_int"] + PRIMARY_COV[coh])
            fj = ols(cc["_y"], Xj)
            e, lo, hi, p = lincom(fj, {"int": 1, "zg_int": -1})
            rows.append(dict(cohort=coh, outcome=lab, model="Joint model contrast", term="pTau217xPVeD minus GFAPxPVeD", n=fj["n"], beta=e, lo=lo, hi=hi, p=p))
            e, lo, hi, p = lincom(fj, {"int": 1, "zn_int": -1})
            rows.append(dict(cohort=coh, outcome=lab, model="Joint model contrast", term="pTau217xPVeD minus NfLxPVeD", n=fj["n"], beta=e, lo=lo, hi=hi, p=p))
    t = pd.DataFrame(rows)
    t["fdr"] = np.nan
    for k, g in t.groupby(["cohort", "model", "term"]):
        t.loc[g.index, "fdr"] = bh(g.p.values)
    t.to_csv(OUT + "T_specificity.csv", index=False)


# ---------------------------------------------------------------- 7. partial volume / atrophy / protocol / stage
def robustness():
    rows = []
    for coh, df in DATA.items():
        lst = [x for x in OUTC[coh] if x[0] in ("Amyloid PET", "Tau PET")]
        specs = [("Primary", PRIMARY_COV[coh], (), None)]
        if coh == "ADNI":
            specs += [("+ CSF fraction (ventricular/atrophy proxy) + ICV", PRIMARY_COV[coh] + ["csf_frac", "zicv"], (), None),
                      ("+ brain parenchymal fraction + hippocampal volume + ICV", PRIMARY_COV[coh] + ["brain_frac", "hippo", "zicv"], (), None),
                      ("+ diffusion QC metrics (QA index, QSDR R)", PRIMARY_COV[coh] + ["qa", "qsdr_R"], (), None),
                      ("+ ADNI phase + DTI protocol", PRIMARY_COV[coh], ("phase", "protocol"), None),
                      ("+ amyloid PET tracer", PRIMARY_COV[coh], ("amy_tracer",), None),
                      ("Cognitively unimpaired only", PRIMARY_COV[coh], (), df.dx == 0),
                      ("MCI or dementia only", PRIMARY_COV[coh], (), df.dx >= 1)]
        else:
            specs += [("+ brain parenchymal fraction (BrainSegNotVent/ICV)", PRIMARY_COV[coh] + ["brain_frac"], (), None),
                      ("+ brain parenchymal fraction + hippocampal volume + ICV", PRIMARY_COV[coh] + ["brain_frac", "hippo", "zicv"], (), None),
                      ("+ diffusion QC metric (QA index)", PRIMARY_COV[coh] + ["qa"], (), None),
                      ("+ DTI protocol", PRIMARY_COV[coh], ("protocol",), None),
                      ("+ ethnicity (Hispanic) and race (Black)", PRIMARY_COV[coh] + ["hispanic", "race_black"], (), None),
                      ("Cognitively unimpaired only", PRIMARY_COV[coh], (), df.dx == 0),
                      ("MCI or dementia only", PRIMARY_COV[coh], (), df.dx >= 1),
                      ("Non-Hispanic participants", PRIMARY_COV[coh], (), df.hispanic == 0),
                      ("Hispanic participants", PRIMARY_COV[coh], (), df.hispanic == 1)]
        for lab, cv, cat, sub in specs:
            for dom, y, olab in lst:
                if lab.startswith("+ amyloid PET tracer") and dom != "Amyloid PET":
                    continue
                _, r = fit_int(df, y, cv, extra_cat=cat, sub=sub)
                rows.append(dict(cohort=coh, analysis=lab, outcome=olab, **r))
        # heterogeneity by clinical stage (int x impaired)
        for dom, y, olab in lst:
            dd = df.copy(); dd["int_imp"] = dd["int"] * dd.impaired; dd["zp_imp"] = dd.zp * dd.impaired; dd["zv_imp"] = dd.zv * dd.impaired
            _, r = fit_int(dd, y, PRIMARY_COV[coh] + ["impaired", "zp_imp", "zv_imp", "int"], term="int_imp")
            rows.append(dict(cohort=coh, analysis="Heterogeneity: interaction x clinical impairment (P for difference)", outcome=olab, **r))
        # PVeD vs atrophy correlations
    t = pd.DataFrame(rows)
    t.to_csv(OUT + "T_robustness.csv", index=False)
    cor = []
    for coh, df in DATA.items():
        for v in ["csf_frac", "brain_frac", "hippo", "wmh", "age", "qa", "qsdr_R"]:
            if v in df and df[v].notna().sum() > 30:
                dd = df[["pved", v]].dropna()
                rr, pp = stats.spearmanr(dd.pved, dd[v])
                cor.append(dict(cohort=coh, variable=v, n=len(dd), spearman_rho=rr, p=pp))
    pd.DataFrame(cor).to_csv(OUT + "T_pved_correlates.csv", index=False)


# ---------------------------------------------------------------- 8. Table 1
def table1():
    order = ["Low pTau217 / High PVeD", "Low pTau217 / Low PVeD", "High pTau217 / High PVeD", "High pTau217 / Low PVeD"]
    rows = []
    for coh, df in DATA.items():
        cont = [("Age, years", "age", "mean"), ("Education, years", "edu", "mean"), ("Plasma pTau217, pg/mL", "ptau_raw", "median"), ("Total PVeD", "pved", "mean"),
                ("CDR sum of boxes", "cdrsb", "median"), ("MMSE", "mmse", "median"), ("PACC", "pacc", "mean"), ("Amyloid PET, Centiloid", "amy_cl", "median"),
                ("Tau PET, meta-temporal SUVR" if coh == "ADNI" else "Tau PET, MsTL tau-stage weighted SUVR", "tau_meta" if coh == "ADNI" else "tau_mstl", "median"),
                ("Plasma GFAP, pg/mL", "gfap_raw", "median"), ("Plasma NfL, pg/mL", "nfl_raw", "median"),
                ("Hippocampal volume, mL", "hippo", "mean"), ("WMH volume (log)", "wmh", "mean")]
        cat = [("Female sex", lambda d: 1 - d.male), ("APOE ε4 carrier", lambda d: d.apoe4), ("Hispanic ethnicity", lambda d: d.hispanic),
               ("Black race", lambda d: d.race_black), ("Amyloid PET positive", lambda d: d.amy_pos)]
        n_row = dict(cohort=coh, variable="N", Total=str(len(df)))
        for o in order:
            n_row[o] = str((df.profile == o).sum())
        rows.append(n_row)
        for lab, v, kind in cont:
            x = df[v].copy()
            if v == "hippo" and coh == "HABS-HD":
                x = x / 1000
            r = dict(cohort=coh, variable=lab, n_avail=int(x.notna().sum()))
            fn = (lambda s: f"{s.mean():.2f} ({s.std():.2f})") if kind == "mean" else (lambda s: f"{s.median():.2f} [{s.quantile(.25):.2f}, {s.quantile(.75):.2f}]")
            if v in ("ptau_raw", "pved", "amy_suvr") or "tau" in v:
                fn = (lambda s: f"{s.mean():.3f} ({s.std():.3f})") if kind == "mean" else (lambda s: f"{s.median():.3f} [{s.quantile(.25):.3f}, {s.quantile(.75):.3f}]")
            r["Total"] = fn(x.dropna())
            groups = []
            for o in order:
                s = x[df.profile == o].dropna(); r[o] = fn(s); groups.append(s)
            p = stats.f_oneway(*groups).pvalue if kind == "mean" else stats.kruskal(*groups).pvalue
            r["P"] = p; r["test"] = "ANOVA" if kind == "mean" else "Kruskal-Wallis"
            rows.append(r)
        for lab, f in cat:
            x = f(df)
            if x.notna().sum() == 0:
                continue
            r = dict(cohort=coh, variable=lab + ", n (%)", n_avail=int(x.notna().sum()))
            r["Total"] = f"{int(x.sum())} ({x.mean() * 100:.1f})"
            tab = []
            for o in order:
                s = x[df.profile == o].dropna(); r[o] = f"{int(s.sum())} ({s.mean() * 100:.1f})"; tab.append([s.sum(), len(s) - s.sum()])
            r["P"] = stats.chi2_contingency(np.array(tab))[1]; r["test"] = "Chi-square"
            rows.append(r)
        # diagnosis
        x = df.dx
        r = dict(cohort=coh, variable="Clinical stage, n (%)", n_avail=int(x.notna().sum()))
        tab = pd.crosstab(df.profile, x)
        r["P"] = stats.chi2_contingency(tab.loc[order].values)[1]; r["test"] = "Chi-square"
        rows.append(r)
        names = {0: "CN/CU", 1: "MCI", 2: "Dementia"}
        for k in [0, 1, 2]:
            rr = dict(cohort=coh, variable="  " + names[k])
            rr["Total"] = f"{int((x == k).sum())} ({(x == k).sum() / x.notna().sum() * 100:.1f})"
            for o in order:
                s = x[df.profile == o]
                rr[o] = f"{int((s == k).sum())} ({(s == k).sum() / s.notna().sum() * 100:.1f})"
            rows.append(rr)
    t = pd.DataFrame(rows)
    t.to_csv(OUT + "T_table1.csv", index=False)
    # missingness of covariates
    miss = []
    for coh, df in DATA.items():
        for v in ["age", "male", "edu", "apoe4", "cdrsb", "dx", "amy_cl", "tau_meta" if coh == "ADNI" else "tau_mstl", "hippo", "wmh", "pacc", "mmse", "gfap_raw", "nfl_raw"]:
            miss.append(dict(cohort=coh, variable=v, available=int(df[v].notna().sum()), missing=int(df[v].isna().sum())))
    pd.DataFrame(miss).to_csv(OUT + "T_missingness.csv", index=False)


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    if step in ("all", "models"):
        models_table(); simple_slopes(); temporal(); print("models done", flush=True)
    if step in ("all", "disc"):
        discrimination(); print("disc done", flush=True)
    if step in ("all", "spec"):
        specificity(); robustness(); table1(); print("spec done", flush=True)
    if step in ("all", "long"):
        longitudinal(); print("long done", flush=True)
