"""ADNI longitudinal PACC: model-form comparison.
Endpoint: the three-test ADNI PACC composite analysed in the original submission (mean of available z-scores for
ADAS-Cog delayed word recall [reversed], MMSE and Logical Memory delayed recall, referenced to ADNI baseline
cognitively normal participants); sensitivity: all three tests required. Rows include
(i) the direct profile-by-time mixed model (revision primary), (ii) the two-stage approach of Huang et al. 2024
(participant slopes from a random intercept/slope LME adjusted for age, sex and education, then regressed on
profiles), (iii) least-squares two-stage slopes, and (iv) reconstruction of the original submission, which used a
three-test composite and empirical-Bayes (BLUP) slopes. Two-stage analyses that treat predicted slopes as data
ignore their uncertainty and understate standard errors."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import pandas as pd, numpy as np, sys, warnings; warnings.filterwarnings('ignore')
from statlib import *; from analyses import A, lmm_three_way
from scipy import stats
RA = CFG.ADNI
L0 = pd.read_csv(RA + "pacc_long_pved.csv"); L0['id'] = L0.RID.astype(str); L0 = L0[L0.id.isin(A.id)]
P = ['Low pTau217 / Low PVeD', 'High pTau217 / High PVeD', 'High pTau217 / Low PVeD']
A['dx_mci'] = (A.dx == 1).astype(float).where(A.dx.notna()); A['dx_dem'] = (A.dx == 2).astype(float).where(A.dx.notna())
cov = ['za', 'male', 'ze']
def long_df(col, filt=None):
    Q = L0 if filt is None else L0[filt]
    Q = Q[['id', 'time_from_dti_years', col]].rename(columns={'time_from_dti_years': 't', col: 'y'}).dropna()
    Q['y'] = Q.y / bsd[col]
    return Q.merge(A, on='id')
bsd = {'PACC_COMPOSITE': A.pacc.std()}
def stage2(S, cv):
    S = S.merge(A, on='id')
    for p in P: S['p_'+p] = (S.profile == p).astype(float)
    f = ols(S.slope, design(S, ['p_'+p for p in P] + cv))
    return f['n'], *lincom(f, {'p_High pTau217 / Low PVeD': 1, 'p_High pTau217 / High PVeD': -1})
def blup_slopes(L, covs=()):
    L = L.dropna(subset=list(covs)); cols = [np.ones(len(L)), L.t.values]
    for c in covs: cols += [L[c].astype(float).values, L[c].astype(float).values * L.t.values]
    X = np.column_stack(cols); Z = X[:, :2]; f = lmm_fit(L.y.values, X, Z, L.id.values)
    D, s2, b = f['D'], f['sigma2'], f['beta']; out = []
    for i, idx in L.groupby('id').indices.items():
        Zi = Z[idx]; Vi = Zi @ D @ Zi.T + s2 * np.eye(len(idx)); r = L.y.values[idx] - X[idx] @ b
        fixed = b[1] + sum(b[3 + 2*k] * L[c].astype(float).values[idx[0]] for k, c in enumerate(covs))
        out.append((i, fixed + (D @ Zi.T @ np.linalg.solve(Vi, r))[1]))
    return pd.DataFrame(out, columns=['id', 'slope'])
def ols_slopes(L):
    return pd.DataFrame([(i, stats.linregress(g.t, g.y).slope) for i, g in L.groupby('id') if len(g) >= 2 and g.t.max() > g.t.min()], columns=['id', 'slope'])
def direct(L, cv):
    L = L.dropna(subset=cv); Xd = {'const': np.ones(len(L)), 't': L.t.values}
    for p in P: z_ = (L.profile == p).astype(float).values; Xd[p] = z_; Xd[p+':t'] = z_ * L.t.values
    for c in cv: Xd[c] = L[c].astype(float).values; Xd[c+':t'] = Xd[c] * L.t.values
    X = pd.DataFrame(Xd); f = lmm_fit(L.y.values, X.values, np.column_stack([np.ones(len(L)), L.t.values]), L.id.values)
    n = list(X.columns); w = np.zeros(len(n)); w[n.index('High pTau217 / Low PVeD:t')] = 1; w[n.index('High pTau217 / High PVeD:t')] = -1
    e = w @ f['beta']; se = np.sqrt(w @ f['cov'] @ w); return L.id.nunique(), e, e-1.96*se, e+1.96*se, 2*stats.norm.sf(abs(e/se))
rows = []
def add(endp, lab, res):
    n, e, lo, hi, p = res; rows.append(dict(cohort='ADNI', endpoint=endp, method=lab, n=n, est=e, lo=lo, hi=hi, p=p))
for endp, col, filt in [('Three-test PACC composite (primary)', 'PACC_COMPOSITE', None),
                        ('Three-test composite, all three tests required', 'PACC_COMPOSITE', L0.N_PACC_COMPONENTS >= 3)]:
    L = long_df(col, filt); c = L.groupby('id').size(); L2 = L[L.id.isin(c[c >= 2].index)]
    add(endp, 'Direct mixed model, profile x time (all participants)', direct(L, cov))
    add(endp, 'Direct mixed model, profile x time (>=2 assessments)', direct(L2, cov))
    r = lmm_three_way(L, cov, endp); add(endp, 'Direct mixed model, continuous time x pTau217 x PVeD', (r['n_participants'], r['three_way'], r['three_way_lo'], r['three_way_hi'], r['three_way_p']))
    if filt is None:
        S = blup_slopes(L, ['age', 'male', 'edu']); S = S[S.id.isin(L2.id)]
        add(endp, 'Two-stage: LME-derived slopes adjusted for age, sex, education (Huang et al. 2024)', stage2(S, cov))
        B = blup_slopes(L); B = B[B.id.isin(L2.id)]
        add(endp, 'Two-stage: empirical-Bayes (BLUP) slopes', stage2(B, cov))
        add(endp, 'Two-stage: empirical-Bayes (BLUP) slopes + APOE e4', stage2(B, cov + ['apoe4']))
        add(endp, 'Two-stage: empirical-Bayes (BLUP) slopes + clinical diagnosis', stage2(B, cov + ['dx_mci', 'dx_dem']))
        add(endp, 'Two-stage: per-participant least-squares slopes', stage2(ols_slopes(L2), cov))
t = pd.DataFrame(rows); t.to_csv(CFG.OUT + 'T_long_method_comparison.csv', index=False); print(t.round(4).to_string())
