#!/usr/bin/env python3
"""Regenerate the five current manuscript figures from audited analysis outputs.

FIG1 estimator validation, empirical SSF, objective-vs-self-report check
FIG2 masking mechanism, recovery fidelity, sparse-sample false positives
FIG3 attenuation identification sensitivity and SSF-calibrated power sensitivity
FIG4 hormone level versus cycle phase
FIG5 instrument/design sensitivity over the unidentified coupled-smooth fraction q
"""
from __future__ import annotations
import argparse,os,warnings
warnings.filterwarnings('ignore')
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h4_frontier as H

def _need(res,names):
    paths=[os.path.join(res,n) for n in names]; missing=[p for p in paths if not os.path.exists(p)]
    if missing: raise FileNotFoundError(', '.join(missing))
    return paths

def fig1(res,out):
    pssf,pobj=_need(res,['table08_instrument_ssf.csv','table14_objective.csv']); ssf=pd.read_csv(pssf); obj=pd.read_csv(pobj); fig,ax=plt.subplots(1,3,figsize=(15.5,4.5)); x=np.arange(3); w=.36
    ax[0].bar(x-w/2,[.077,.036,.028],w,label='mean |bias|'); ax[0].bar(x+w/2,[.157,.082,.082],w,label='max |bias|'); ax[0].set_xticks(x); ax[0].set_xticklabels(['AR(1)','ACF-linear','Spectral']); ax[0].set_ylabel('bias'); ax[0].set_title('(a) Estimator validation'); ax[0].legend(fontsize=8)
    vals=ssf.spectral.values; ax[1].bar(np.arange(len(vals)),vals); ax[1].set_xticks(np.arange(len(vals))); ax[1].set_xticklabels([str(m)[:14] for m in ssf.measure],rotation=40,ha='right',fontsize=7); ax[1].set_ylabel('smooth-signal fraction'); ax[1].set_title('(b) Empirical SSF')
    xx=np.arange(len(obj)); ax[2].bar(xx-.18,obj.SD_ri,.36,label='observed SD(r_i)'); ax[2].bar(xx+.18,obj.null_SD,.36,label='surrogate-null median')
    for i,p in enumerate(obj.p): ax[2].text(i,max(obj.SD_ri.iloc[i],obj.null_SD.iloc[i])+.004,f'p={p:.3f}',ha='center',fontsize=7)
    ax[2].set_xticks(xx); ax[2].set_xticklabels(obj.outcome,rotation=20,ha='right',fontsize=8); ax[2].set_title('(c) Objective vs self-report (exploratory)'); ax[2].legend(fontsize=7); plt.tight_layout(); path=f'{out}/FIG1_estimators_ssf_objective.png'; plt.savefig(path,dpi=220,bbox_inches='tight'); plt.close(); return path

def fig2(res,out):
    p3,p4=_need(res,['table03_null_calibration.csv','table04_fidelity.csv']); nul=pd.read_csv(p3); fid=pd.read_csv(p4); fig,ax=plt.subplots(1,3,figsize=(15.5,4.5)); xg=np.linspace(.15,.95,300); ax[0].plot(xg,H.inverted_u(xg),lw=2); mu=H.DA_OPT-H.K_GAIN*H._E2_MEAN
    for b in (mu-.16,mu+.16):
        lo,hi=b+H.K_GAIN*.13,b+H.K_GAIN; xs=np.linspace(lo,hi,60); ax[0].plot(xs,H.inverted_u(xs),lw=5)
    ax[0].axvline(H.DA_OPT,ls=':'); ax[0].set_xlabel('dopaminergic tone'); ax[0].set_ylabel('performance'); ax[0].set_title('(a) Bidirectional masking mechanism')
    piv=fid.pivot(index='obs_per_person',columns='sigma_b',values='fidelity').sort_index(); im=ax[1].imshow(piv.values,origin='lower',aspect='auto',vmin=0,vmax=1); ax[1].set_xticks(range(len(piv.columns))); ax[1].set_xticklabels([f'{v:.2f}' for v in piv.columns]); ax[1].set_yticks(range(len(piv.index))); ax[1].set_yticklabels(piv.index); ax[1].set_xlabel('heterogeneity (sigma_b)'); ax[1].set_ylabel('observations/person'); ax[1].set_title('(b) Recovery fidelity'); plt.colorbar(im,ax=ax[1],label='corr(estimated,true)')
    ax[2].plot(nul.obs_per_person,nul.crit1_fp,'o-',label='effect-size criterion alone'); ax[2].plot(nul.obs_per_person,nul.lrt_fp,'o-',label='LRT alone'); ax[2].plot(nul.obs_per_person,nul.both_fp,'o-',label='dual rule'); ax[2].axhline(.05,ls='--',label='alpha=.05'); ax[2].set_ylim(-.03,1.05); ax[2].set_xlabel('observations/person'); ax[2].set_ylabel('false-positive rate'); ax[2].set_title('(c) Why sparse thresholds fail'); ax[2].legend(fontsize=7); plt.tight_layout(); path=f'{out}/FIG2_mechanism_fidelity.png'; plt.savefig(path,dpi=220,bbox_inches='tight'); plt.close(); return path

def fig3(res,out):
    p9,pp=_need(res,['table09_attenuation_sensitivity.csv','registered_test_power_sensitivity.csv']); att=pd.read_csv(p9); pw=pd.read_csv(pp); fig,ax=plt.subplots(1,3,figsize=(15.5,4.4)); ax[0].plot(att.f,att.attenuation,'o-'); ax[0].set_xlabel('f: high-frequency variance treated as genuine'); ax[0].set_ylabel('implied attenuation'); ax[0].set_ylim(0,1.05); ax[0].set_title('(a) Identification sensitivity')
    for q,g in pw.groupby('coupled_fraction'): ax[1].plot(g.sigma_b,g.power,'o-',label=f'q={q:.2f}')
    ax[1].axhline(.80,ls='--'); ax[1].set_xlabel('sigma_b'); ax[1].set_ylabel('surrogate-test power'); ax[1].set_ylim(-.03,1.03); ax[1].set_title('(b) Power not identified by SSF alone'); ax[1].legend(fontsize=7)
    piv=pw.pivot(index='coupled_fraction',columns='sigma_b',values='power').sort_index(); im=ax[2].imshow(piv.values,origin='lower',aspect='auto',vmin=0,vmax=1); ax[2].set_xticks(range(len(piv.columns))); ax[2].set_xticklabels([f'{v:.3g}' for v in piv.columns],rotation=30); ax[2].set_yticks(range(len(piv.index))); ax[2].set_yticklabels([f'{v:.2f}' for v in piv.index]); ax[2].set_xlabel('sigma_b'); ax[2].set_ylabel('q'); ax[2].set_title('(c) SSF-calibrated sensitivity surface'); plt.colorbar(im,ax=ax[2],label='power'); plt.tight_layout(); path=f'{out}/FIG3_attenuation_power_sensitivity.png'; plt.savefig(path,dpi=220,bbox_inches='tight'); plt.close(); return path

def fig4(res,out):
    (p,)=_need(res,['fig04_phase_locked.csv']); d=pd.read_csv(p); fig,ax=plt.subplots(1,2,figsize=(10.5,4.3)); x=np.arange(len(d)); w=.38; ax[0].bar(x-w/2,np.abs(d.r_e2_level),w,label='|r| with E2 level'); ax[0].bar(x+w/2,d.phase_eta2,w,label='eta^2 of cycle phase'); ax[0].set_xticks(x); ax[0].set_xticklabels(d.item,rotation=25,ha='right'); ax[0].legend(fontsize=8); ax[0].set_title('(a) Level vs phase'); ax[1].bar(x,d.menstrual_mean); ax[1].axhline(0,lw=1); ax[1].set_xticks(x); ax[1].set_xticklabels(d.item,rotation=25,ha='right'); ax[1].set_ylabel('menstrual-phase within-person mean'); ax[1].set_title('(b) Menstrual-phase shift'); plt.tight_layout(); path=f'{out}/FIG4_phase_locked.png'; plt.savefig(path,dpi=220,bbox_inches='tight'); plt.close(); return path

def fig5(res,out):
    pg,pb=_need(res,['h4v2_ssf_sensitivity_grid.csv','budget_sensitivity.csv']); grid=pd.read_csv(pg); bud=pd.read_csv(pb); fig,ax=plt.subplots(1,3,figsize=(16,4.5)); uq=sorted(grid.coupled_fraction.unique()); q0=float(uq[min(1,len(uq)-1)]); sb0=float(sorted(grid.sigma_b.unique())[-1]); sub=grid[(grid.coupled_fraction==q0)&(grid.sigma_b==sb0)]; piv=sub.pivot(index='ssf_outcome',columns='ssf_predictor',values='dual_rate').sort_index(); im=ax[0].imshow(piv.values,origin='lower',aspect='auto',vmin=0,vmax=1); ax[0].set_xticks(range(len(piv.columns))); ax[0].set_xticklabels([f'{v:.2f}' for v in piv.columns]); ax[0].set_yticks(range(len(piv.index))); ax[0].set_yticklabels([f'{v:.2f}' for v in piv.index]); ax[0].set_xlabel('predictor SSF'); ax[0].set_ylabel('outcome SSF'); ax[0].set_title(f'(a) Instrument sensitivity (q={q0:.2f})'); plt.colorbar(im,ax=ax[0],label='dual rate')
    scen=list(dict.fromkeys(bud.scenario.tolist())); qvals=sorted(bud.coupled_fraction.unique())
    for q in qvals:
        g=bud[bud.coupled_fraction==q].set_index('scenario').reindex(scen); ax[1].plot(range(len(scen)),g.dual_rate,'o-',label=f'q={q:.2f}'); ax[2].plot(range(len(scen)),g.fidelity,'o-',label=f'q={q:.2f}')
    labels=[s.replace(' observations','') for s in scen]; ax[1].set_xticks(range(len(scen))); ax[1].set_xticklabels(labels,rotation=45,ha='right',fontsize=7); ax[1].set_ylim(-.03,1.03); ax[1].set_ylabel('dual-criterion rate'); ax[1].set_title('(b) Budget sensitivity'); ax[1].legend(fontsize=6); ax[2].axhline(.70,ls='--'); ax[2].set_xticks(range(len(scen))); ax[2].set_xticklabels(labels,rotation=45,ha='right',fontsize=7); ax[2].set_ylim(0,1.03); ax[2].set_ylabel('recovery fidelity'); ax[2].set_title('(c) Person-specific recovery'); plt.tight_layout(); path=f'{out}/FIG5_design_sensitivity.png'; plt.savefig(path,dpi=220,bbox_inches='tight'); plt.close(); return path

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--results',default='./results'); ap.add_argument('--out',default='./figures'); a=ap.parse_args(); os.makedirs(a.out,exist_ok=True); made=[]
    for name,fn in [('FIG1',fig1),('FIG2',fig2),('FIG3',fig3),('FIG4',fig4),('FIG5',fig5)]:
        try: made.append(fn(a.results,a.out)); print(f'[made {name}] {made[-1]}')
        except FileNotFoundError as e: print(f'[skip {name}] missing prerequisite: {e}')
    print(f'[done] {len(made)}/5 figures regenerated')
if __name__=='__main__': main()
