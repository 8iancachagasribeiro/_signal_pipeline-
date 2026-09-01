from __future__ import annotations
import argparse, os, time
import numpy as np, pandas as pd
import h4_frontier as H
from brm_common import *

NS=(10,20,30,40,56); CYCLES=(1,2,3,4)
SCENARIOS=[
 ('reference','reference',{}),
 ('nonstationarity','linear_drift',{'drift':'linear'}),
 ('nonstationarity','slow_baseline',{'drift':'slow'}),
 ('residual_AR1','rho_0.2',{'resid_rho':.2}),
 ('residual_AR1','rho_0.5',{'resid_rho':.5}),
 ('residual_AR1','rho_0.8',{'resid_rho':.8}),
 ('heteroscedasticity','luteal_1.5x',{'hetero':True}),
 ('functional_form','sigmoid',{'form':'sigmoid'}),
 ('functional_form','asym_linear',{'form':'asym_linear'}),
 ('functional_form','piecewise',{'form':'piecewise'}),
 ('individual_gain','K_sd_0.05',{'vary_K':True}),
 ('directional_distribution','shifted_unimodal',{'direction':'shifted'}),
 ('directional_distribution','bimodal',{'direction':'bimodal'}),
 ('measurement_error','white_rel_0.55',{'reliability':.55}),
 ('measurement_error','AR1_rho_0.3',{'meas_rho':.3}),
 ('contamination','outliers_2pct',{'outliers':.02}),
 ('missingness','MCAR_10pct',{'missing':'mcar10'}),
 ('missingness','MCAR_25pct',{'missing':'mcar25'}),
 ('missingness','phase_dependent',{'missing':'phase'}),
 ('missingness','outcome_MNAR_like',{'missing':'mnar'}),
]

def simulate_cell(seed,n_obs,n_cycles,opts,R=1000,N=H.N_SUBJ):
    rng=np.random.default_rng(seed); span=H.CYCLE_LEN*n_cycles
    days=np.linspace(0,span,n_obs,endpoint=False); x=H.e2(days)
    form=opts.get('form','inverted_u')
    direction=opts.get('direction','balanced')
    if direction=='shifted': b=rng.normal(BASE_B+.05,.10,size=(R,N))
    elif direction=='bimodal':
        comp=rng.random((R,N))<.5; b=BASE_B+np.where(comp,-.11,.11)+rng.normal(0,.035,size=(R,N))
    else: b=rng.normal(BASE_B,.10,size=(R,N))
    if opts.get('vary_K'):
        K=np.clip(rng.normal(H.K_GAIN,.05,size=(R,N)),.02,.30)
    else: K=np.full((R,N),H.K_GAIN)
    sig=response(b[...,None]+K[...,None]*x[None,None,:],form)
    drift=opts.get('drift')
    if drift=='linear':
        sig=sig + np.linspace(-.075,.075,n_obs)[None,None,:]
    elif drift=='slow':
        sig=sig + .055*np.sin(np.linspace(0,1.25*np.pi,n_obs))[None,None,:]
    rho=opts.get('resid_rho',0.0)
    if opts.get('hetero'):
        mult=np.where((days%28)>=15,1.5,1.0); mult=mult/np.sqrt(np.mean(mult*mult))
        eps=rng.normal(size=sig.shape)*H.SIGMA_STATE*mult[None,None,:]
    else: eps=ar1_noise(rng,sig.shape,rho,H.SIGMA_STATE)
    y_true=sig+eps
    rel=opts.get('reliability',REL_REF); v=np.var(y_true,axis=(1,2),keepdims=True)
    sd=np.sqrt(np.maximum(v,1e-12)*(1-rel)/rel)
    meas_rho=opts.get('meas_rho',0.0)
    meas=ar1_noise(rng,sig.shape,meas_rho,1.0)*sd
    y=y_true+meas
    if opts.get('outliers'):
        p=float(opts['outliers']); mask=rng.random(y.shape)<p
        scale=np.std(y,axis=(1,2),keepdims=True)
        y=y+mask*rng.normal(0,5.0,size=y.shape)*scale
    miss=opts.get('missing')
    if miss:
        if miss=='mcar10': pm=np.full(y.shape,.10)
        elif miss=='mcar25': pm=np.full(y.shape,.25)
        elif miss=='phase':
            base=np.where((days%28)>=15,.30,.05); pm=np.broadcast_to(base,y.shape)
        elif miss=='mnar':
            zy=(y-y.mean(axis=-1,keepdims=True))/np.maximum(y.std(axis=-1,keepdims=True),1e-9)
            pm=.04+.24/(1+np.exp(-1.3*zy))
        mask=rng.random(y.shape)<pm; y=np.where(mask,np.nan,y)
    est=corr_rows(y,np.broadcast_to(x,y.shape))
    truth=true_coupling(b,K if opts.get('vary_K') else None,form)
    mets=rep_metrics(est,truth); out=summarize_metrics(mets,R)
    out.update(n_obs=n_obs,n_cycles=n_cycles,total_span_days=span)
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--reps',type=int,default=1000); ap.add_argument('--out-dir',default='./results'); ap.add_argument('--seed',type=int,default=SEEDS['robustness']); ap.add_argument('--scenario-start',type=int,default=0); ap.add_argument('--scenario-end',type=int,default=None); ap.add_argument('--tag',default='')
    a=ap.parse_args(); os.makedirs(a.out_dir,exist_ok=True); rows=[]; k=0; t0=time.time()
    end=len(SCENARIOS) if a.scenario_end is None else a.scenario_end
    selected=list(enumerate(SCENARIOS))[a.scenario_start:end]
    for gi,(group,level,opts) in selected:
      for cyc in CYCLES:
       for n in NS:
        seed=a.seed+gi*10000+cyc*100+n; rec=simulate_cell(seed,n,cyc,opts,a.reps); rec.update(factor=group,level=level,seed=seed); rows.append(rec); k+=1
        if k%20==0: print(f'[robustness] {k}/{len(selected)*len(CYCLES)*len(NS)} cells',flush=True)
    df=pd.DataFrame(rows); path=os.path.join(a.out_dir,f'brm_robustness_metrics{a.tag}.csv'); df.to_csv(path,index=False)
    print(f'[saved] {path}; cells={len(df)}; elapsed={time.time()-t0:.1f}s')
if __name__=='__main__': main()
