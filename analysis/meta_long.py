"""Two-cohort random-effects (DerSimonian-Laird) synthesis of longitudinal PACC estimates (baseline-SD units/year)."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import pandas as pd, numpy as np
from scipy import stats
c = pd.read_csv(CFG.OUT + 'T_lmm_continuous.csv'); p = pd.read_csv(CFG.OUT + 'T_lmm_profiles.csv')
def dl(est, se):
    est, se = np.asarray(est), np.asarray(se); w = 1/se**2; fe = (w*est).sum()/w.sum()
    Q = (w*(est-fe)**2).sum(); k = len(est); tau2 = max(0, (Q-(k-1))/(w.sum()-(w**2).sum()/w.sum()))
    wr = 1/(se**2+tau2); re = (wr*est).sum()/wr.sum(); sre = np.sqrt(1/wr.sum())
    I2 = max(0, (Q-(k-1))/Q)*100 if Q > 0 else 0
    return re, re-1.96*sre, re+1.96*sre, 2*stats.norm.sf(abs(re/sre)), fe, Q, stats.chi2.sf(Q, k-1), I2, tau2
rows = []
for lab, mc, mp in [('age, sex, education', 'LMM-1: age, sex, education', 'age, sex, education'), ('+ APOE e4', 'LMM-2: + APOE e4', '+ APOE e4'), ('+ APOE e4 + CDR-SB', 'LMM-3: + APOE e4 + baseline CDR-SB', '+ APOE e4 + CDR-SB')]:
    g = c[c.model == mc].set_index('cohort')
    for estimand, est, lo, hi in [('Continuous time x pTau217 x PVeD', 'three_way', 'three_way_lo', 'three_way_hi')]:
        e = g[est].values; se = (g[hi]-g[lo]).values/3.92
        r = dl(e, se); rows.append(dict(estimand=estimand, model=lab, ADNI=f"{e[0]:.3f}", HABS=f"{e[1]:.3f}", pooled=r[0], lo=r[1], hi=r[2], p=r[3], Q_p=r[6], I2=r[7], tau2=r[8]))
    h = p[(p.model == mp) & (p.term.str.startswith('Contrast'))].set_index('cohort')
    e = h.est.values; se = (h.hi-h.lo).values/3.92; r = dl(e, se)
    rows.append(dict(estimand='Profile contrast (high/low vs high/high PVeD)', model=lab, ADNI=f"{e[0]:.3f}", HABS=f"{e[1]:.3f}", pooled=r[0], lo=r[1], hi=r[2], p=r[3], Q_p=r[6], I2=r[7], tau2=r[8]))
t = pd.DataFrame(rows); t.to_csv(CFG.OUT + 'T_long_meta.csv', index=False); print(t.round(4).to_string())
