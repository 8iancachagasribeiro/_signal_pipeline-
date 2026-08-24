#!/usr/bin/env python3
"""Design and instrument-quality sensitivity with predictor and outcome SSF.

Empirical SSF is not injected as psychometric reliability. Predictor and outcome series
are calibrated to requested spectral SSF. Outcome SSF does not identify how much smooth
variance is attributable to E2, so ``coupled_fraction`` q is explicit and results are
sensitivity analyses rather than point estimates of empirical power.
"""
from __future__ import annotations
import argparse, os, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import h4_frontier as H
from run_h4 import within_r_matrix
from fastlrt import lrt_random_slope_varying_x
from ssf_estimators import ssf_spectral
from ssf_power import add_white_noise_for_estimated_ssf,smooth_nuisance,_standardize
N_SUBJ=42; SEED=17
_XD=H.e2(np.linspace(0,H.CYCLE_LEN,400,endpoint=False))

def true_coupling(b): return within_r_matrix(H.inverted_u(b[:,None]+H.K_GAIN*_XD[None,:]),_XD)

def simulate(rng,n_obs,n_cycles,ssf_x,ssf_y,sigma_b,coupled_fraction,n_subj=N_SUBJ):
    span=int(np.ceil(H.CYCLE_LEN*n_cycles))
    days=np.linspace(0,H.CYCLE_LEN*n_cycles,n_obs,endpoint=False) if n_obs>span else np.sort(rng.choice(np.arange(span),size=n_obs,replace=False)).astype(float)
    x_true=H.e2(days); b=rng.normal(H.DA_OPT-H.K_GAIN*H._E2_MEAN,sigma_b,n_subj)
    Y=np.empty((n_subj,len(days))); X=np.empty_like(Y)
    for i in range(n_subj):
        mech=H.inverted_u(b[i]+H.K_GAIN*x_true); m=_standardize(mech); z=smooth_nuisance(len(days),rng)
        smooth_y=_standardize(np.sqrt(coupled_fraction)*m+np.sqrt(1-coupled_fraction)*z)
        Y[i]=add_white_noise_for_estimated_ssf(smooth_y,ssf_y,rng,ssf_spectral)
        X[i]=add_white_noise_for_estimated_ssf(x_true,ssf_x,rng,ssf_spectral)
    return Y,X,b

def eval_cell(rng,n_obs,n_cycles,ssf_x,ssf_y,sigma_b,coupled_fraction,n_sims=200):
    dual=lrt=0; fid=[]; medr=[]
    for _ in range(n_sims):
        yo,xo,b=simulate(rng,n_obs,n_cycles,ssf_x,ssf_y,sigma_b,coupled_fraction)
        rh=np.array([np.corrcoef(xo[i],yo[i])[0,1] for i in range(len(yo))]); rt=true_coupling(b)
        prop=float(np.nanmean(np.abs(rh)>.20)); p=lrt_random_slope_varying_x(yo,xo)
        if np.isfinite(p):
            if p<.05:
                lrt+=1
                if prop>.50: dual+=1
        m=np.isfinite(rh)&np.isfinite(rt)
        if m.sum()>3 and np.std(rh[m])>1e-9 and np.std(rt[m])>1e-9: fid.append(np.corrcoef(rh[m],rt[m])[0,1])
        medr.append(np.nanmedian(np.abs(rh)))
    return dual/n_sims,lrt/n_sims,float(np.mean(fid)) if fid else np.nan,float(np.mean(medr))

def instrument_grid(rng,n_sims,q_grid):
    rows=[]
    for q in q_grid:
        for sb in [.10,.20]:
            for sx in [.35,.469,.60,.75,.90]:
                for sy in [.25,.323,.50,.70,.90]:
                    d,l,f,mr=eval_cell(rng,85,3,sx,sy,sb,q,n_sims)
                    rows.append(dict(coupled_fraction=q,sigma_b=sb,ssf_predictor=sx,ssf_outcome=sy,n_obs=85,n_cycles=3,dual_rate=d,lrt_rate=l,fidelity=f,median_abs_r=mr))
                    print(f'q={q:.2f} sb={sb:.2f} SSFx={sx:.3f} SSFy={sy:.3f} dual={d:.2f} fid={f:.2f}',flush=True)
    return pd.DataFrame(rows)

def budget_sensitivity(rng,n_sims,q_grid):
    px,py=.469,.323
    scenarios=[('baseline 85 obs / 3 cycles',px,py,85,3),('2x observations / 6 cycles',px,py,170,6),('4x observations / 12 cycles',px,py,340,12),('improve outcome SSF to .70',px,.70,85,3),('improve predictor SSF to .70',.70,py,85,3),('improve both SSF to .70',.70,.70,85,3),('improve both SSF to .90',.90,.90,85,3)]
    rows=[]
    for q in q_grid:
        for name,sx,sy,nobs,ncyc in scenarios:
            d,l,f,mr=eval_cell(rng,nobs,ncyc,sx,sy,.15,q,n_sims)
            rows.append(dict(coupled_fraction=q,scenario=name,ssf_predictor=sx,ssf_outcome=sy,n_obs=nobs,n_cycles=ncyc,dual_rate=d,lrt_rate=l,fidelity=f,median_abs_r=mr))
            print(f'q={q:.2f} {name}: dual={d:.2f}, fidelity={f:.2f}')
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('mode',nargs='?',default='grid',choices=['grid','budget']); ap.add_argument('--n-sims',type=int,default=200); ap.add_argument('--q',default='0.10,0.25,0.50,0.75,1.00'); ap.add_argument('--seed',type=int,default=SEED); ap.add_argument('--out-dir',default='./results')
    a=ap.parse_args(); os.makedirs(a.out_dir,exist_ok=True); rng=np.random.default_rng(a.seed); qs=[float(v) for v in a.q.split(',')]
    if a.mode=='grid': instrument_grid(rng,a.n_sims,qs).to_csv(f'{a.out_dir}/h4v2_ssf_sensitivity_grid.csv',index=False)
    else: budget_sensitivity(rng,a.n_sims,qs).to_csv(f'{a.out_dir}/budget_sensitivity.csv',index=False)
    print('\nThese are sensitivity analyses over q, not point estimates of empirical power.')
if __name__=='__main__': main()
