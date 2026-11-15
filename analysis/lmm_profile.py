import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import pandas as pd,numpy as np,sys,warnings; warnings.filterwarnings('ignore')
from statlib import *; from analyses import load_long, A, H
from scipy import stats
LA,LH=load_long()
def prof_lmm(L,base_sd,covs,lab,coh):
    L=L.copy(); L['y']=(L.y-0)/base_sd
    L=L.dropna(subset=['y','t']+covs)
    P=['Low pTau217 / Low PVeD','High pTau217 / High PVeD','High pTau217 / Low PVeD']
    Xd={'const':np.ones(len(L)),'t':L.t.values}
    for p in P:
        Xd[p]=(L.profile==p).astype(float).values; Xd[p+':t']=Xd[p]*L.t.values
    for c in covs: Xd[c]=L[c].astype(float).values; Xd[c+':t']=Xd[c]*L.t.values
    X=pd.DataFrame(Xd); f=lmm_fit(L.y.values,X.values,np.column_stack([np.ones(len(L)),L.t.values]),L.id.values)
    names=list(X.columns); b=pd.Series(f['beta'],names); C=pd.DataFrame(f['cov'],names,names)
    def lc(w):
        v=np.array([w.get(k,0) for k in names]); e=v@b.values; se=np.sqrt(v@C.values@v); return e,e-1.96*se,e+1.96*se,2*stats.norm.sf(abs(e/se))
    rows=[]
    # covariates centred? za,ze standardized, male not centred -> report slope at male=mean
    mm={c+':t':L[c].mean() for c in covs}
    rows.append(dict(cohort=coh,model=lab,term='Low pTau217 / High PVeD',**dict(zip(['est','lo','hi','p'],lc({'t':1,**mm})))))
    for p in P: rows.append(dict(cohort=coh,model=lab,term=p,**dict(zip(['est','lo','hi','p'],lc({'t':1,p+':t':1,**mm})))))
    rows.append(dict(cohort=coh,model=lab,term='Contrast: High/Low vs High/High',**dict(zip(['est','lo','hi','p'],lc({'High pTau217 / Low PVeD:t':1,'High pTau217 / High PVeD:t':-1})))))
    for r in rows: r.update(n=L.id.nunique(),obs=len(L))
    return rows
rows=[]
sdA=1.0; sdH=1.0  # outcomes already in baseline-SD units (load_long)
for coh,L,sd in [('ADNI',LA,sdA),('HABS-HD',LH,sdH)]:
    for lab,cv,dat in [('age, sex, education',['za','male','ze'],L),('+ APOE e4',['za','male','ze','apoe4'],L),('+ APOE e4 + CDR-SB',['za','male','ze','apoe4','zcdr'],L),('>=2 assessments, age/sex/education',['za','male','ze'],L[L.id.isin(L.groupby('id').size().loc[lambda s:s>=2].index)])]:
        if coh=='HABS-HD': dat=dat[dat.t>=-0.5]
        rows+=prof_lmm(dat,sd,cv,lab,coh); print(coh,lab,rows[-1]['est'],rows[-1]['p'],flush=True)
pd.DataFrame(rows).to_csv(CFG.OUT + 'T_lmm_profiles.csv',index=False)
print('SD',sdA,sdH)
