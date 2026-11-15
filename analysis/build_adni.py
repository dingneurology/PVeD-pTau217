import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import pandas as pd,numpy as np,sys,warnings; warnings.filterwarnings("ignore")
from statlib import *
R=CFG.ADNI
b=pd.read_csv(R+'ADNI_PVED_AandD_github_package_20260702/results/alz20_pacc_reanalysis/baseline_main_with_alz20_pacc.csv',low_memory=False)
b['DTI_DATE']=pd.to_datetime(b.DTI_DATE)
def nearest(src,cols,name,date='SCANDATE'):
    src=src.copy(); src[date]=pd.to_datetime(src[date],errors='coerce')
    m=b[['RID','DTI_DATE']].merge(src[['RID',date]+cols],on='RID')
    m['dd']=(m[date]-m.DTI_DATE).dt.days; m['ad']=m.dd.abs()
    m=m.dropna(subset=['dd']).sort_values('ad').drop_duplicates('RID')
    return m.drop(columns=['DTI_DATE',date]).rename(columns={c:name+'_'+c for c in cols+['dd','ad']})
a=pd.read_csv(R+'UCBERKELEY_AMY_6MM_28Apr2026.csv',low_memory=False)
t=pd.read_csv(R+'UCBERKELEY_TAU_6MM_28Apr2026.csv',low_memory=False)
w=pd.read_csv(R+'UCD_WMH_24May2026.csv',low_memory=False)
c=pd.read_csv(R+'CDR_28Apr2026.csv',low_memory=False)
tcols=['META_TEMPORAL_SUVR','CTX_ENTORHINAL_SUVR','CTX_INFERIORTEMPORAL_SUVR','CTX_MIDDLETEMPORAL_SUVR','CTX_FUSIFORM_SUVR','CTX_PRECUNEUS_SUVR','CTX_POSTERIORCINGULATE_SUVR']
tcols=[x for x in tcols if x in t.columns]
A=nearest(a,['CENTILOIDS','SUMMARY_SUVR','AMYLOID_STATUS','TRACER'],'amy')
T=nearest(t,tcols+['TRACER'],'tau')
W=nearest(w,['TOTAL_HIPPO','TOTAL_WMH','CEREBRUM_TCV','TOTAL_BRAIN','TOTAL_CSF'],'mri','EXAMDATE')
C=nearest(c,['CDGLOBAL','CDRSB'],'cdr','VISDATE')
d=b.merge(A,on='RID',how='left').merge(T,on='RID',how='left').merge(W,on='RID',how='left').merge(C,on='RID',how='left')
for col in ['cdr_CDGLOBAL','cdr_CDRSB']:
    d.loc[d[col]<0,col]=np.nan
d['ptau']=np.log10(d.bio_plasma_ptau217.where(d.bio_plasma_ptau217>0))
d['z_ptau']=z(d.ptau); d['z_pved']=z(d.PVeD_total); d['z_age']=z(d.AGE); d['z_edu']=z(d.EDUCATION)
d['male']=(d.SEX==1).astype(float)
d['apoe4']=d.APOE4_carrier
d['int']=d.z_ptau*d.z_pved
d['hippo_icv']=d.mri_TOTAL_HIPPO/d.mri_CEREBRUM_TCV
d['log_wmh']=np.log(d.mri_TOTAL_WMH+1)
d['log_wmh_icv']=np.log(d.mri_TOTAL_WMH/d.mri_CEREBRUM_TCV*100+0.01)
for m in ['gfap','nfl']:
    v=d['bio_plasma_'+m]; d['z_'+m]=z(np.log10(v.where(v>0)))
d.to_pickle(CFG.WORK + 'adni.pkl')
print(d.shape)
