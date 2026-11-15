import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import pandas as pd,numpy as np,sys,warnings; warnings.filterwarnings('ignore')
from statlib import *
R=CFG.HABS
def nid(s): return s.astype(str).str.extract(r'^(\d+)',expand=False).fillna(s.astype(str)).str.replace(r'\.0$','',regex=True)
def num(s): 
    x=pd.to_numeric(s,errors='coerce'); return x.where(~x.isin([-9999,-8888,-777777,-888888,-999999]))
def pdate(s):
    n=pd.to_numeric(s,errors='coerce'); ex=pd.to_datetime('1899-12-30')+pd.to_timedelta(n,unit='D'); tx=pd.to_datetime(s,errors='coerce'); return tx.mask(n.notna(),ex)
def mdate(y,m): return pd.to_datetime(dict(year=pd.to_numeric(y,errors='coerce'),month=pd.to_numeric(m,errors='coerce'),day=15),errors='coerce')
def incl(f,kw='Inclusion'):
    xl=pd.ExcelFile(R+f); parts=[]
    for s in xl.sheet_names:
        if kw.lower() in s.lower() and 'fail' not in s.lower() and 'exclu' not in s.lower():
            d=pd.read_excel(R+f,sheet_name=s)
            if 'Med_ID' in d and len(d): d['sheet']=s; d['Med_ID']=nid(d.Med_ID); parts.append(d)
    return pd.concat(parts,ignore_index=True)
pv=pd.read_csv(R+'PVeD_metrics.csv'); pv['Med_ID']=nid(pv.Med_ID); pv=pv.drop_duplicates('Med_ID')
mr=pd.read_csv(R+'HABS_HD_ADNI_mirror_tables_figures_20260703/validation/HABS_pacc_like_cross_sectional_model_ready.csv',low_memory=False)
mr['Med_ID']=nid(mr.Med_ID); mr=mr.drop_duplicates('Med_ID')[['Med_ID','pved_mri_date']]; mr['pved_mri_date']=pd.to_datetime(mr.pved_mri_date)
cl=pd.read_csv(R+'RP_HD_7_Clinical.csv',low_memory=False); cl['Med_ID']=nid(cl.Med_ID); cl['Visit_ID']=num(cl.Visit_ID); cl['vdate']=pdate(cl.Visit_Date)
bio=pd.read_excel(R+'RP_HD_7_Biomarkers.xlsx',sheet_name=0); bio['Med_ID']=nid(bio.Med_ID); bio['Visit_ID']=num(bio.Visit_ID)
bio=bio.merge(cl[['Med_ID','Visit_ID','vdate']],on=['Med_ID','Visit_ID'],how='left')
ge=pd.read_excel(R+'RP_HD_7_Genomics.xlsx'); ge['Med_ID']=nid(ge.Med_ID); ge=ge.sort_values(['Med_ID','Visit_ID']).drop_duplicates('Med_ID')
d=pv.merge(mr,on='Med_ID',how='left')
def firstbio(col,name):
    b=bio[['Med_ID','Visit_ID','vdate',col]].copy(); b[col]=num(b[col]); b=b[b[col].notna()].sort_values(['Med_ID','vdate','Visit_ID']).drop_duplicates('Med_ID')
    return b[['Med_ID',col,'vdate','Visit_ID']].rename(columns={col:name,'vdate':name+'_date','Visit_ID':name+'_visit'})
for col,name in [('r7_QTX_Plasma_pTau217_PLUS','ptau_plus'),('r7_QTX_Plasma_pTau217_V2','ptau_v2'),('r7_QTX_Plasma_GFAP','gfap_r7'),('r7_QTX_Plasma_GFAP_PLUS','gfap_plus'),('r7_QTX_Plasma_NfL_V2','nfl_r7'),('r5_QTX_Plasma_GFAP','gfap_r5'),('r5_QTX_Plasma_NfL','nfl_r5'),('r3_QTX_Plasma_NfL','nfl_r3')]:
    d=d.merge(firstbio(col,name),on='Med_ID',how='left')
# baseline clinical (first visit) + clinical at plasma visit
clin_cols=['Med_ID','Visit_ID','vdate','Age','ID_Gender','ID_Education','CDX_Cog','CDR_Global','CDR_Sum','MMSE_Total','ID_Race','ID_Hispanic']
clin_cols=[c for c in clin_cols if c in cl.columns]
cb=cl.sort_values(['Med_ID','Visit_ID']).drop_duplicates('Med_ID')[clin_cols]
d=d.merge(cb.rename(columns={'Visit_ID':'clin_visit','vdate':'clin_date'}),on='Med_ID',how='left')
# clinical at the plasma visit (for CDR)
cp=cl[['Med_ID','Visit_ID','CDR_Global','CDR_Sum','CDX_Cog','Age']].rename(columns={'Visit_ID':'ptau_plus_visit','CDR_Global':'CDR_Global_pv','CDR_Sum':'CDR_Sum_pv','CDX_Cog':'CDX_Cog_pv','Age':'Age_pv'})
d=d.merge(cp.drop_duplicates(['Med_ID','ptau_plus_visit']),on=['Med_ID','ptau_plus_visit'],how='left')
d=d.merge(ge[['Med_ID','APOE4_Positivity','APOE4_Genotype','ID_Education']].rename(columns={'ID_Education':'edu_ge'}),on='Med_ID',how='left')
amy=incl('HABS_HD_RP_8_FBB_PET_data_BraskieLab_20260501.xlsx'); amy['amy_date']=mdate(amy.FBB_PET_year,amy.FBB_PET_month)
for c in ['Global_FBB_SUVR','Centiloid_Units_GlobalAmyloid','AB_pos_FBB']: amy[c]=num(amy[c])
tau=incl('HABS_HD_RP_8_PI2620_PET_Data_BraskieLab_20260501.xlsx'); tau['tau_date']=mdate(tau.PI2620_PET_year,tau.PI2620_PET_month)
tc=['Medial_Temporal_PI2620_SUVR','Posterior_Cingulate_PI2620_SUVR','Lateral_Parietal_PI2620_SUVR','Tennant_MsTL_taucutpoint_weighted_SUVR']
for c in tc: tau[c]=num(tau[c])
hip=incl('HABS_HD_RP_8_Hippocampal_Volume_HippoDeep_Data_BraskieLab_20260501.xlsx','Inclusions'); hip['hip_date']=mdate(hip.MRI_year,hip.MRI_month)
for c in ['L_hippocampus','R_hippocampus','ICV_HD']: hip[c]=num(hip[c])
ct=incl('HABS_HD_RP_8_CorticalThickness_ATN_MetaROI_Data_BraskieLab_20260501.xlsx','Inclusions'); ct['ct_date']=mdate(ct.MRI_year,ct.MRI_month); ct['CorticalThickness_MetaROI']=num(ct.CorticalThickness_MetaROI)
xl=pd.ExcelFile(R+'HABS_HD_RP_8_WMH_Combination_Report_Analysis_Borzage_Lab_20260507.xlsx')
wm=pd.concat([pd.read_excel(xl,sheet_name=s) for s in xl.sheet_names if 'fail' not in s]); wm['Med_ID']=nid(wm.Med_ID); wm['wmh_date']=mdate(wm.MRI_Date_Year,wm.MRI_Date_Month); wm['wmh_log']=num(wm['02_WMH_Volume_Log'])
def pick(src,dcol,cols,pref,mode):
    m=d[['Med_ID','pved_mri_date']].merge(src[['Med_ID',dcol]+cols],on='Med_ID')
    m['mo']=(m[dcol]-m.pved_mri_date).dt.days/30.4375; m['amo']=m.mo.abs()
    m=m.dropna(subset=[cols[0]])
    if mode=='first': m=m.sort_values(['Med_ID',dcol]).drop_duplicates('Med_ID')
    else: m=m.sort_values(['Med_ID','amo']).drop_duplicates('Med_ID')
    return m[['Med_ID',dcol]+cols+['mo']].rename(columns={**{c:pref+c for c in cols},dcol:pref+'date','mo':pref+'mo'})
for mode in ['first','near']:
    d=d.merge(pick(amy,'amy_date',['Centiloid_Units_GlobalAmyloid','Global_FBB_SUVR','AB_pos_FBB'],f'A{mode}_',mode),on='Med_ID',how='left')
    d=d.merge(pick(tau,'tau_date',tc,f'T{mode}_',mode),on='Med_ID',how='left')
d=d.merge(pick(hip,'hip_date',['L_hippocampus','R_hippocampus','ICV_HD'],'H_','near'),on='Med_ID',how='left')
d=d.merge(pick(ct,'ct_date',['CorticalThickness_MetaROI'],'CT_','near'),on='Med_ID',how='left')
d=d.merge(pick(wm,'wmh_date',['wmh_log'],'W_','near'),on='Med_ID',how='left')
d['ptau_mo']=(d.ptau_plus_date-d.pved_mri_date).dt.days/30.4375
d.to_pickle(CFG.WORK + 'habs_raw.pkl')
# DTI acquisition protocol per scan (used for the acquisition-protocol sensitivity analysis)
pd.read_excel(R + 'HABS_HD_RP_8_DTI_ROI_Results.xlsx', usecols=['Med_ID', 'Year', 'Month', 'Protocol']).to_pickle(CFG.WORK + 'habs_dti.pkl'); print(d.shape)
print(d[['PVeD_total','ptau_plus','ptau_v2','gfap_r7','gfap_plus','nfl_r7','gfap_r5','nfl_r5','nfl_r3','Afirst_Centiloid_Units_GlobalAmyloid','Tfirst_Medial_Temporal_PI2620_SUVR','H_L_hippocampus','W_wmh_log','CT_CorticalThickness_MetaROI','CDR_Global','CDR_Sum_pv','pved_mri_date']].notna().sum())
