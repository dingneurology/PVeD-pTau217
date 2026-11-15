import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import pandas as pd,numpy as np,sys,warnings; warnings.filterwarnings('ignore')
from statlib import *; from analyses import OUTC,DATA,PRIMARY_COV
rows=[]
for coh,lst in OUTC.items():
    df=DATA[coh]
    for dom,y,lab in lst:
        if dom not in ('Amyloid PET','Tau PET'): continue
        d=df.copy(); d['_y']=z(d[y]); cc=d.dropna(subset=['_y','zp','zv']+PRIMARY_COV[coh])
        fits={k:ols(cc._y,design(cc,v)) for k,v in {'cov':PRIMARY_COV[coh],'p':PRIMARY_COV[coh]+['zp'],'pv':PRIMARY_COV[coh]+['zp','zv'],'int':PRIMARY_COV[coh]+['zp','zv','int']}.items()}
        rows.append(dict(cohort=coh,outcome=lab,n=fits['p']['n'],adjR2_cov=fits['cov']['adj_r2'],adjR2_ptau=fits['p']['adj_r2'],adjR2_ptau_pved=fits['pv']['adj_r2'],adjR2_int=fits['int']['adj_r2'],
            d_int_vs_ptau=fits['int']['adj_r2']-fits['p']['adj_r2'],aic_ptau=fits['p']['aic'],aic_int=fits['int']['aic'],bic_ptau=fits['p']['bic'],bic_int=fits['int']['bic']))
t=pd.DataFrame(rows); t.to_csv(CFG.OUT + 'T_incremental_R2.csv',index=False); print(t.round(3).to_string())
