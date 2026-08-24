#!/usr/bin/env python3
"""Null calibration, recovery fidelity, and aliasing for the theoretical H4 design."""
import argparse, warnings
import numpy as np, pandas as pd
import h4_frontier as H
from fastlrt import lrt_random_slope_fast
warnings.filterwarnings('ignore')
N_SUBJ=H.N_SUBJ; REL_REF=.73
_XD=H.e2(np.linspace(0,H.CYCLE_LEN,400,endpoint=False))

def within_r(Y,x):
    xc=x-x.mean(); Yc=Y-Y.mean(axis=1,keepdims=True); den=np.sqrt((Yc**2).sum(axis=1)*(xc**2).sum()); den=np.where(den<1e-12,np.nan,den); return (Yc@xc)/den

def true_coupling(b): return within_r(H.inverted_u(b[:,None]+H.K_GAIN*_XD[None,:]),_XD)

def simulate(rng,days,reliability,sigma_b,n_subj=N_SUBJ):
    x=H.e2(np.asarray(days,float)); b=rng.normal(H.DA_OPT-H.K_GAIN*H._E2_MEAN,sigma_b,n_subj); signal=H.inverted_u(b[:,None]+H.K_GAIN*x[None,:]); y_true=signal+rng.normal(0,H.SIGMA_STATE,signal.shape)
    sd=np.sqrt(max(y_true.var(),1e-12)*(1-reliability)/reliability); return y_true+rng.normal(0,sd,y_true.shape),x,b

def even_days(obs_per_cycle,n_cycles): return np.concatenate([np.linspace(0,H.CYCLE_LEN,obs_per_cycle,endpoint=False)+c*H.CYCLE_LEN for c in range(n_cycles)])
PHASE_TARGETED={2:[3,13],3:[3,13,21],4:[3,9,13,21],5:[3,9,13,17,21],7:[2,6,9,13,17,21,25],10:list(np.linspace(1,27,10))}

def null_calibration(rng,n_sims=500):
    rows=[]
    for d in (2,3,5,7,10,14,21,28):
        c1=c2=both=valid=0
        for _ in range(n_sims):
            Y,x,_=simulate(rng,even_days(d,2),REL_REF,.001); prop=float(np.nanmean(np.abs(within_r(Y,x))>.20)); p=lrt_random_slope_fast(Y,x)
            if not np.isfinite(p): continue
            valid+=1; a=prop>.50; b=p<.05; c1+=a; c2+=b; both+=(a and b)
        den=max(valid,1); rows.append(dict(obs_per_person=d*2,crit1_fp=c1/den,lrt_fp=c2/den,both_fp=both/den,n_valid=valid))
    out=pd.DataFrame(rows); print(out.to_string(index=False)); print('\nEffect-size thresholds alone fail at very sparse sampling; no universal four-observation detection threshold is claimed.'); return out

def fidelity(rng,n_sims=500):
    rows=[]
    for d in (2,3,5,7,10,14,21,28):
        for sb in (.05,.10,.15,.20):
            vals=[]
            for _ in range(n_sims):
                Y,x,b=simulate(rng,even_days(d,2),REL_REF,sb); rh,rt=within_r(Y,x),true_coupling(b); m=np.isfinite(rh)&np.isfinite(rt)
                if m.sum()>3 and np.std(rh[m])>1e-9 and np.std(rt[m])>1e-9: vals.append(np.corrcoef(rh[m],rt[m])[0,1])
            rows.append(dict(obs_per_cycle=d,obs_per_person=d*2,sigma_b=sb,fidelity=float(np.mean(vals)) if vals else np.nan))
    out=pd.DataFrame(rows); print(out.pivot(index='obs_per_person',columns='sigma_b',values='fidelity').round(2)); return out

def aliasing(rng,n_sims=500,sigma_b=.10):
    rows=[]
    for k in (2,3,4,5,7,10):
        out={}
        for tag,days1 in (('even',even_days(k,1)),('phase',np.asarray(PHASE_TARGETED[k],float))):
            days=np.concatenate([days1+c*H.CYCLE_LEN for c in range(2)]); vals=[]
            for _ in range(n_sims):
                Y,x,b=simulate(rng,days,REL_REF,sigma_b); rh,rt=within_r(Y,x),true_coupling(b); m=np.isfinite(rh)&np.isfinite(rt)
                if m.sum()>3 and np.std(rh[m])>1e-9: vals.append(np.corrcoef(rh[m],rt[m])[0,1])
            out[tag]=(float(np.mean(vals)) if vals else np.nan,float(H.e2(days).max()-H.e2(days).min()))
        rows.append(dict(obs_per_cycle=k,fidelity_even=out['even'][0],fidelity_phase=out['phase'][0],e2_range_even=out['even'][1],e2_range_phase=out['phase'][1]))
    out=pd.DataFrame(rows); print(out.to_string(index=False)); return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out-dir',default='./results'); ap.add_argument('--seed',type=int,default=99); ap.add_argument('--n-sims',type=int,default=500,help='Monte Carlo replicates per cell')
    a=ap.parse_args(); import os; os.makedirs(a.out_dir,exist_ok=True); rng=np.random.default_rng(a.seed)
    null_calibration(rng,a.n_sims).to_csv(f'{a.out_dir}/table03_null_calibration.csv',index=False); fidelity(rng,a.n_sims).to_csv(f'{a.out_dir}/table04_fidelity.csv',index=False); aliasing(rng,a.n_sims).to_csv(f'{a.out_dir}/table05_aliasing.csv',index=False); print(f'[saved] {a.out_dir}/')
if __name__=='__main__': main()
