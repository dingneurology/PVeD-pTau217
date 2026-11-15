"""Export source data for the revised Figure 3 and Supplementary Figure S1 (plotted in R: rebuild_revised_Figure3_FigS1.R)."""
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import pandas as pd, numpy as np, sys, warnings; warnings.filterwarnings('ignore')
from statlib import *; from analyses import DATA, PRIMARY_COV, load_long
O = CFG.OUT; S = CFG.FIG_SOURCE
short = {'Amyloid PET Centiloid': 'Amyloid Centiloid', 'Amyloid PET SUVR': 'Amyloid SUVR (within-tracer z)',
         'Tau PET meta-temporal SUVR': 'Tau meta-temporal', 'Tau PET entorhinal SUVR': 'Tau entorhinal', 'Tau PET inferior temporal SUVR': 'Tau inferior temporal',
         'Tau PET middle temporal SUVR': 'Tau middle temporal', 'Tau PET fusiform SUVR': 'Tau fusiform', 'Tau PET precuneus SUVR': 'Tau precuneus',
         'Tau PET posterior cingulate SUVR': 'Tau posterior cingulate', 'Hippocampal volume': 'Hippocampal volume', 'WMH volume (log)': 'WMH',
         'Baseline PACC': 'Baseline PACC', 'Baseline MMSE': 'Baseline MMSE', 'FBB amyloid Centiloid': 'Amyloid Centiloid', 'FBB global amyloid SUVR': 'Amyloid SUVR',
         'PI-2620 medial temporal SUVR': 'Tau medial temporal', 'PI-2620 posterior cingulate SUVR': 'Tau posterior cingulate', 'PI-2620 lateral parietal SUVR': 'Tau lateral parietal',
         'PI-2620 MsTL tau-stage weighted SUVR': 'Tau MsTL stage-weighted', 'AD-signature cortical thickness': 'AD-signature thickness'}
m = pd.read_csv(O + 'T_models_M1_M4.csv'); m = m[m.model.str.startswith('M1')].copy()
m['outcome_short'] = m.outcome.map(short); m['domain_group'] = np.where(m.domain.isin(['Amyloid PET', 'Tau PET']), 'AD PET pathology', 'Other outcomes')
m['sig'] = m.fdr < 0.05; m['order'] = m.groupby('cohort').cumcount()
m[['cohort', 'domain', 'outcome', 'outcome_short', 'domain_group', 'n', 'beta', 'lo', 'hi', 'p', 'fdr', 'sig', 'order']].to_csv(S + 'Fig3_forest.csv', index=False)
# brain maps (Desikan-Killiany regions; FDR-significant negative interactions only)
regmap = {'ADNI': {'Tau PET entorhinal SUVR': ['entorhinal'], 'Tau PET inferior temporal SUVR': ['inferior temporal'], 'Tau PET middle temporal SUVR': ['middle temporal'],
                   'Tau PET fusiform SUVR': ['fusiform'], 'Tau PET precuneus SUVR': ['precuneus'], 'Tau PET posterior cingulate SUVR': ['posterior cingulate']},
          'HABS-HD': {'PI-2620 medial temporal SUVR': ['entorhinal', 'parahippocampal'], 'PI-2620 posterior cingulate SUVR': ['posterior cingulate'],
                      'PI-2620 lateral parietal SUVR': ['inferior parietal', 'superior parietal'], 'PI-2620 MsTL tau-stage weighted SUVR': ['entorhinal', 'fusiform', 'inferior temporal', 'middle temporal']}}
rows = []
for coh, mp in regmap.items():
    g = m[m.cohort == coh].set_index('outcome')
    for oc, regs in mp.items():
        r = g.loc[oc]
        for rg in regs: rows.append(dict(cohort=coh, region=rg, outcome=oc, beta=r.beta, fdr=r.fdr, sig=bool(r.sig and r.beta < 0)))
b = pd.DataFrame(rows)
b = b.groupby(['cohort', 'region']).apply(lambda x: pd.Series(dict(brain_beta=x.beta.min(), any_fdr_sig=bool(x.sig.any()), outcomes='; '.join(x.outcome)))).reset_index()
b.to_csv(S + 'Fig3_brainmap.csv', index=False)
# model-predicted standardized outcomes at PVeD -1/+1 SD
pred = []
for coh, y, lab in [('ADNI', 'amy_cl', 'Amyloid PET Centiloid'), ('ADNI', 'tau_meta', 'Tau PET meta-temporal SUVR (flortaucipir)'),
                    ('HABS-HD', 'amy_cl', 'FBB amyloid Centiloid'), ('HABS-HD', 'tau_mstl', 'PI-2620 MsTL stage-weighted tau SUVR')]:
    df = DATA[coh].copy(); df['_y'] = z(df[y]); cv = PRIMARY_COV[coh]
    X = design(df, ['zp', 'zv', 'int'] + cv); f = ols(df['_y'], X); b_ = f['coef'].beta; V = f['cov']
    d = df.loc[f['idx']]; grid = np.linspace(d.zp.quantile(.05), d.zp.quantile(.95), 120)
    for zv in (-1, 1):
        for zp_ in grid:
            x = pd.Series(0.0, index=b_.index); x['const'] = 1; x['zp'] = zp_; x['zv'] = zv; x['int'] = zp_ * zv
            for c in cv: x[c] = d[c].astype(float).mean()
            e = float(x @ b_); se = float(np.sqrt(x @ V.values @ x)); tc = stats.t.ppf(.975, f['dfr'])
            pred.append(dict(cohort=coh, outcome=lab, z_ptau217=zp_, pved_level='Lower PVeD (-1 SD)' if zv < 0 else 'Higher PVeD (+1 SD)', fit=e, lwr=e - tc * se, upr=e + tc * se))
pd.DataFrame(pred).to_csv(S + 'Fig3_prediction.csv', index=False)
# sensitivity: temporal alignment and adjustment
tw = pd.read_csv(O + 'T_time_windows.csv'); mm = pd.read_csv(O + 'T_models_M1_M4.csv')
keep = {'ADNI': ['Amyloid PET Centiloid', 'Tau PET meta-temporal SUVR', 'Tau PET fusiform SUVR'], 'HABS-HD': ['FBB amyloid Centiloid', 'PI-2620 posterior cingulate SUVR', 'PI-2620 MsTL tau-stage weighted SUVR']}
lab = {'All available': 'Primary (all available)', '±24 months': '±24 months', '±12 months': '±12 months', '±6 months': '±6 months', 'All available + interval covariates': '+ interval covariates'}
rows = []
for coh, outs in keep.items():
    for _, r in tw[(tw.cohort == coh) & tw.outcome.isin(outs)].iterrows():
        rows.append(dict(cohort=coh, outcome=r.outcome, outcome_short=short[r.outcome], model=lab[r.window], n=r.n, beta=r.beta, lo=r.lo, hi=r.hi, fdr=r.fdr))
    for mod, ml in [('M2', 'M2 harmonized covariates'), ('M3', 'M3 + CDR-SB')]:
        for _, r in mm[(mm.cohort == coh) & mm.model.str.startswith(mod) & mm.outcome.isin(outs)].iterrows():
            rows.append(dict(cohort=coh, outcome=r.outcome, outcome_short=short[r.outcome], model=ml, n=r.n, beta=r.beta, lo=r.lo, hi=r.hi, fdr=r.fdr))
pd.DataFrame(rows).to_csv(S + 'Fig3_sensitivity.csv', index=False)
# Supplementary Figure S1
LA, LH = load_long()
sc, pc, ep = [], [], []
for coh, df, L in [('ADNI', DATA['ADNI'], LA), ('HABS-HD', DATA['HABS-HD'], LH)]:
    for _, r in df.iterrows(): sc.append(dict(cohort=coh, log10_ptau217=r.ptau_log, pved_total=r.pved, profile=r.profile))
    for p_, n_ in df.profile.value_counts().items(): pc.append(dict(cohort=coh, profile=p_, n=n_))
    Lp = L if coh == 'ADNI' else L[L.t >= -0.5]
    two = int((Lp.groupby('id').size() >= 2).sum())
    if coh == 'ADNI':
        lst = [('PVeD + plasma\npTau217', len(df)), ('Amyloid PET', df.amy_cl.notna().sum()), ('Meta-temporal tau\n(flortaucipir)', df.tau_meta.notna().sum()), ('Baseline PACC', df.pacc.notna().sum()), ('Longitudinal\nPACC (≥2)', two)]
    else:
        lst = [('PVeD + plasma\npTau217', len(df)), ('Amyloid PET', df.amy_cl.notna().sum()), ('MsTL tau', df.tau_mstl.notna().sum()), ('Baseline PACC', df.pacc.notna().sum()), ('Longitudinal\nPACC (≥2)', two)]
    for k, (l, n_) in enumerate(lst): ep.append(dict(cohort=coh, order=k, label=l, n=int(n_)))
    med = dict(cohort=coh, ptau_median=df.ptau_log.median(), pved_median=df.pved.median())
    pc.append(dict(cohort=coh, profile='__medians__', n=np.nan, **{k: v for k, v in med.items() if k != 'cohort'}))
pd.DataFrame(sc).to_csv(S + 'FigS1_scatter.csv', index=False); pd.DataFrame(pc).to_csv(S + 'FigS1_profile_counts.csv', index=False); pd.DataFrame(ep).to_csv(S + 'FigS1_endpoints.csv', index=False)
print(pd.DataFrame(ep).to_string()); print(pd.DataFrame(pc).to_string()); print(b.to_string())

# ---------------- Figure 2 (profile phenotypes) and Figure 1 replacement data panels
rows = []
dxlab = {'ADNI': {0: 'CN', 1: 'MCI', 2: 'AD dementia'}, 'HABS-HD': {0: 'CU', 1: 'MCI', 2: 'Dementia'}}
for coh, df in DATA.items():
    tau = df.tau_meta if coh == 'ADNI' else df.tau_mstl
    for i, r in df.iterrows():
        rows.append(dict(cohort=coh, profile=r.profile, dx=dxlab[coh].get(r.dx, np.nan) if pd.notna(r.dx) else np.nan,
                         amyloid_centiloid=r.amy_cl, tau_suvr=tau.loc[i], log10_gfap=np.log10(r.gfap_raw) if r.gfap_raw > 0 else np.nan,
                         log10_nfl=np.log10(r.nfl_raw) if r.nfl_raw > 0 else np.nan, pacc=r.pacc))
pd.DataFrame(rows).to_csv(S + 'Fig2_participants.csv', index=False)
p = pd.read_csv(O + 'T_lmm_profiles.csv'); p = p[(p.model == 'age, sex, education')].copy()
p.to_csv(S + 'Fig1_profile_slopes.csv', index=False)
print('Fig2/Fig1 sources written')
