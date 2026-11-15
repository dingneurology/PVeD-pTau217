import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_ROOT = _os.path.dirname(_HERE)
_sys.path[:0] = [_HERE, _ROOT, _os.path.join(_ROOT, "analysis")]
import config as CFG  # noqa: E402
import pandas as pd,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#555','axes.linewidth':0.6})
O=CFG.OUT; F=CFG.FIG
col={'Low pTau217 / High PVeD':'#2C6FBB','Low pTau217 / Low PVeD':'#8C8C8C','High pTau217 / High PVeD':'#E08A1E','High pTau217 / Low PVeD':'#C0392B'}
order=list(col)
# ---- Figure 4 revised: LMM profile slopes
t=pd.read_csv(O+'T_lmm_profiles.csv'); t=t[t.model=='age, sex, education']
fig,ax=plt.subplots(1,2,figsize=(7.2,3.0))
for a,coh,lab in [(ax[0],'ADNI','A  ADNI (PACC, baseline-SD units/year)'),(ax[1],'HABS-HD','B  HABS-HD (PACC, baseline-SD units/year)')]:
    s=t[t.cohort==coh].set_index('term')
    for i,p in enumerate(order):
        r=s.loc[p]; a.errorbar(i,r.est,yerr=[[r.est-r.lo],[r.hi-r.est]],fmt='o',color=col[p],ms=6,capsize=0,lw=2)
    c=s.loc['Contrast: High/Low vs High/High']
    a.axhline(0,color='#999',lw=0.6,ls='--')
    a.set_xticks(range(4)); a.set_xticklabels([p.replace(' / ','\n') for p in order],fontsize=7)
    a.set_title(lab,loc='left',fontsize=8.5,fontweight='bold')
    a.text(0.02,0.04,f"High/Low vs High/High PVeD:\n{c.est:.3f} ({c.lo:.3f} to {c.hi:.3f}); P = {c.p:.3f}",transform=a.transAxes,fontsize=7)
    a.set_ylabel('Adjusted annual change (95% CI)')
fig.tight_layout(); fig.savefig(F+'Figure_4_revised_LMM_profiles.pdf'); fig.savefig(F+'Figure_4_revised_LMM_profiles.png',dpi=300)
# ---- Supplementary Figure S2: simple slopes
s=pd.read_csv(O+'T_simple_slopes.csv'); s=s[~s.outcome.str.contains('SUVR') | s.outcome.str.contains('Tau|PI-2620')]
s=s[s.outcome!='FBB global amyloid SUVR']
fig,ax=plt.subplots(1,2,figsize=(7.2,3.8),sharex=True)
for a,coh in zip(ax,['ADNI','HABS-HD']):
    d=s[s.cohort==coh].reset_index(drop=True)
    for i,r in d.iterrows():
        for k,off,c,m in [('low',-0.15,'#C0392B','PVeD −1 SD'),('high',0.15,'#2C6FBB','PVeD +1 SD')]:
            a.errorbar(r[f'slope_{k}'],i+off,xerr=[[r[f'slope_{k}']-r[f'slope_{k}_lo']],[r[f'slope_{k}_hi']-r[f'slope_{k}']]],fmt='o',color=c,ms=4,lw=1.5,label=m if i==0 else None)
    a.set_yticks(range(len(d))); a.set_yticklabels([o.replace('Tau PET ','Tau ').replace('PI-2620 ','Tau ').replace(' SUVR','') for o in d.outcome],fontsize=7); a.invert_yaxis()
    a.set_title(('A  ' if coh=='ADNI' else 'B  ')+coh,loc='left',fontweight='bold'); a.set_xlabel('pTau217–PET slope (SD per SD, 95% CI)')
    a.axvline(0,color='#999',lw=0.6,ls='--')
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=2,frameon=False,fontsize=7)
fig.tight_layout(rect=(0,0.06,1,1)); fig.savefig(F+'Supplementary_Figure_S2_simple_slopes.pdf'); fig.savefig(F+'Supplementary_Figure_S2_simple_slopes.png',dpi=300)
# ---- Supplementary Figure S3: time windows
w=pd.read_csv(O+'T_time_windows.csv'); w=w[w.outcome.isin(['Amyloid PET Centiloid','Tau PET meta-temporal SUVR','Tau PET fusiform SUVR','FBB amyloid Centiloid','PI-2620 posterior cingulate SUVR','PI-2620 MsTL tau-stage weighted SUVR'])]
wins=['All available','±24 months','±12 months','±6 months','All available + interval covariates']; wc=dict(zip(wins,['#222','#6A51A3','#2C6FBB','#1B9E77','#999']))
fig,ax=plt.subplots(1,2,figsize=(7.2,3.6),sharex=True)
for a,coh in zip(ax,['ADNI','HABS-HD']):
    d=w[w.cohort==coh]; outs=list(dict.fromkeys(d.outcome)); ws=[x for x in wins if x in set(d.window)]
    for i,o in enumerate(outs):
        for j,x in enumerate(ws):
            r=d[(d.outcome==o)&(d.window==x)].iloc[0]; y=i+(j-(len(ws)-1)/2)*0.15
            a.errorbar(r.beta,y,xerr=[[r.beta-r.lo],[r.hi-r.beta]],fmt='o',ms=3.5,lw=1.2,color=wc[x],label=f'{x}' if i==0 else None)
    a.set_yticks(range(len(outs))); a.set_yticklabels([o.replace('Tau PET ','Tau ').replace('PI-2620 ','Tau ').replace(' SUVR','') for o in outs],fontsize=7); a.invert_yaxis()
    a.axvline(0,color='#999',lw=0.6,ls='--'); a.set_xlabel('pTau217 × PVeD interaction β (95% CI)'); a.set_title(('A  ' if coh=='ADNI' else 'B  ')+coh,loc='left',fontweight='bold')
h,l=ax[1].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=5,frameon=False,fontsize=7)
fig.tight_layout(rect=(0,0.08,1,1)); fig.savefig(F+'Supplementary_Figure_S3_time_windows.pdf'); fig.savefig(F+'Supplementary_Figure_S3_time_windows.png',dpi=300)
print('ok')
