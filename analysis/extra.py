import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import pandas as pd,numpy as np,sys,warnings; warnings.filterwarnings('ignore')
from statlib import *; from analyses import OUTC,DATA,PRIMARY_COV,fit_int
rows=[]
for coh,lst in OUTC.items():
    df=DATA[coh]
    for dom,y,lab in lst:
        if dom not in ('Amyloid PET','Tau PET'): continue
        _,r=fit_int(df,y,PRIMARY_COV[coh],sub=(df.amy_pos==1)); rows.append(dict(cohort=coh,analysis='Amyloid-PET-positive participants',outcome=lab,**r))
        for side,v in [('Left PVeD','pved_L'),('Right PVeD','pved_R')]:
            d=df.copy(); d['zv']=z(d[v]); d['int']=d.zp*d.zv; _,r=fit_int(d,y,PRIMARY_COV[coh]); rows.append(dict(cohort=coh,analysis=side+' as modifier',outcome=lab,**r))
# HABS V2 assay
H=DATA['HABS-HD'].copy(); raw=pd.read_pickle(CFG.WORK + 'habs_raw.pkl')[['Med_ID','ptau_v2']].rename(columns={'Med_ID':'id'})
H=H.merge(raw,on='id',how='left'); H['zp']=z(np.log10(H.ptau_v2.where(H.ptau_v2>0))); H['int']=H.zp*H.zv
for dom,y,lab in OUTC['HABS-HD']:
    if dom in ('Amyloid PET','Tau PET'):
        _,r=fit_int(H,y,PRIMARY_COV['HABS-HD']); rows.append(dict(cohort='HABS-HD',analysis='Alternative Simoa pTau217 (V2) assay',outcome=lab,**r))
Hh=DATA['HABS-HD']
for v,lab in [('pacc','PACC primary (>=3/4 components)'),('pacc_strict','PACC strict (4/4 components)'),('pacc_trails','PACC with Trails B instead of animal fluency'),('pacc_np','PACC neuropsychological components only')]:
    _,r=fit_int(Hh,v,PRIMARY_COV['HABS-HD']); rows.append(dict(cohort='HABS-HD',analysis='Baseline PACC definitions',outcome=lab,**r))
Aa=DATA['ADNI']
for v,lab in [('pacc','Three-test PACC composite (primary; >=2 of 3 tests in 99.7%)'),('pacc_3comp','Three-test composite, all three tests required')]:
    _,r=fit_int(Aa,v,PRIMARY_COV['ADNI']); rows.append(dict(cohort='ADNI',analysis='Baseline PACC definitions',outcome=lab,**r))
t=pd.DataFrame(rows); t['fdr']=np.nan
for k,g in t.groupby(['cohort','analysis']): t.loc[g.index,'fdr']=bh(g.p.values)
t.to_csv(CFG.OUT + 'T_extra_sensitivity.csv',index=False); print(t[['cohort','analysis','outcome','n','beta','lo','hi','p','fdr']].round(3).to_string())
