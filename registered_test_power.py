#!/usr/bin/env python3
"""Power/sensitivity analysis for the preregistered phase-randomized surrogate test.

This file deliberately does NOT report a single "power with the real instruments".
Empirical SSF values (E3G ~= .469; confirmatory self-report ~= .323) identify how much
variance is smooth, but they do not identify what share of the outcome's smooth variance
is actually driven by E3G. Power is therefore evaluated over an explicit sensitivity
parameter ``q`` (coupled_fraction), the share of smooth outcome variance attributable to
the mechanistic E2 -> performance component.

The prior pipeline incorrectly injected SSF values as if they were classical
reliabilities. That implementation is retained nowhere in the inferential path.
"""
from __future__ import annotations
import argparse, os, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import h4_frontier as H
from ssf_estimators import ssf_spectral
from ssf_power import add_white_noise_for_estimated_ssf, smooth_nuisance, summarize_actual_ssf, _standardize

N_SUBJ=42; N_OBS=85; SPAN=90; SSF_X=.469; SSF_Y=.323; B_SURR=500; N_REP=200; ALPHA=.05


def phase_randomize(x,rng):
    x=np.asarray(x,float); n=len(x); X=np.fft.rfft(x); mag=np.abs(X)
    ph=rng.uniform(0,2*np.pi,len(X)); ph[0]=0.0
    if n%2==0: ph[-1]=0.0
    return np.fft.irfft(mag*np.exp(1j*ph),n)


def _corr(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if np.std(a)<1e-12 or np.std(b)<1e-12: return 0.0
    return float(np.corrcoef(a,b)[0,1])


def simulate_study(rng,sigma_b,ssf_x=SSF_X,ssf_y=SSF_Y,coupled_fraction=.5,n_subj=N_SUBJ,n_obs=N_OBS,span=SPAN):
    """Simulate on a regular calendar, then mask to the observed paired days."""
    b=rng.normal(H.DA_OPT-H.K_GAIN*H._E2_MEAN,sigma_b,n_subj)
    offs=rng.uniform(0,H.CYCLE_LEN,n_subj); records=[]
    for i in range(n_subj):
        days_full=np.arange(span,dtype=float); x_true=H.e2(days_full+offs[i])
        mech=H.inverted_u(b[i]+H.K_GAIN*x_true)
        x_obs_full=add_white_noise_for_estimated_ssf(x_true,ssf_x,rng,ssf_spectral)
        q=float(coupled_fraction); m=_standardize(mech); z=smooth_nuisance(len(m),rng)
        smooth_y=_standardize(np.sqrt(q)*m+np.sqrt(1-q)*z)
        y_obs_full=add_white_noise_for_estimated_ssf(smooth_y,ssf_y,rng,ssf_spectral)
        obs_idx=np.sort(rng.choice(np.arange(span),size=n_obs,replace=False))
        records.append({"x_full":x_obs_full,"y_full":y_obs_full,"obs_idx":obs_idx,
                        "x":x_obs_full[obs_idx],"y":y_obs_full[obs_idx]})
    return records


def surrogate_test(records,rng,B=B_SURR):
    r_obs=np.array([_corr(r['x'],r['y']) for r in records]); S_obs=float(np.std(r_obs)); S_null=np.empty(B)
    for k in range(B):
        rs=[]
        for r in records:
            xr_full=phase_randomize(r['x_full'],rng); rs.append(_corr(xr_full[r['obs_idx']],r['y']))
        S_null[k]=np.std(rs)
    p=(1+int(np.sum(S_null>=S_obs)))/(B+1)
    return p,S_obs,float(np.median(S_null))


def ssf_diagnostic(rng,q=.5):
    rec=simulate_study(rng,.10,coupled_fraction=q)
    return summarize_actual_ssf([r['x_full'] for r in rec],ssf_spectral), summarize_actual_ssf([r['y_full'] for r in rec],ssf_spectral)


def run_calibration(rng,reps,surr,q):
    rej=0
    for _ in range(reps):
        rec=simulate_study(rng,1e-6,coupled_fraction=q); p,_,_=surrogate_test(rec,rng,B=surr); rej+=p<ALPHA
    return rej/reps


def run_power(rng,reps,surr,q_grid,sigma_grid):
    rows=[]
    for q in q_grid:
        for sb in sigma_grid:
            rej=0; s_obs=[]
            for _ in range(reps):
                rec=simulate_study(rng,sb,coupled_fraction=q); p,s,_=surrogate_test(rec,rng,B=surr)
                rej+=p<ALPHA; s_obs.append(s)
            rows.append(dict(coupled_fraction=q,sigma_b=sb,power=rej/reps,mean_SD_ri=float(np.mean(s_obs)),reps=reps,surrogates=surr,target_ssf_predictor=SSF_X,target_ssf_outcome=SSF_Y))
            print(f"q={q:.2f} sigma_b={sb:.3f} -> power={rej/reps:.3f}",flush=True)
    return pd.DataFrame(rows)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('mode',nargs='?',default='power',choices=['power','calib','diagnostic'])
    ap.add_argument('--reps',type=int,default=N_REP); ap.add_argument('--surr',type=int,default=B_SURR); ap.add_argument('--seed',type=int,default=23)
    ap.add_argument('--q',default='0.10,0.25,0.50,0.75,1.00'); ap.add_argument('--sigma',default='0.05,0.075,0.10,0.15,0.20,0.30'); ap.add_argument('--out-dir',default='./results')
    a=ap.parse_args(); os.makedirs(a.out_dir,exist_ok=True); rng=np.random.default_rng(a.seed); t0=time.time(); qs=[float(x) for x in a.q.split(',')]
    if a.mode=='diagnostic':
        for q in qs:
            sx,sy=ssf_diagnostic(rng,q=q); print(f"q={q:.2f}: target SSF x/y={SSF_X:.3f}/{SSF_Y:.3f}; median estimated SSF={sx:.3f}/{sy:.3f}")
        return
    if a.mode=='calib':
        rows=[]
        for q in qs:
            fpr=run_calibration(rng,a.reps,a.surr,q); print(f"q={q:.2f}: false-positive rate={fpr:.3f} (alpha={ALPHA:.2f})")
            rows.append(dict(coupled_fraction=q,false_positive_rate=fpr,reps=a.reps,surrogates=a.surr))
        pd.DataFrame(rows).to_csv(f"{a.out_dir}/registered_test_calibration.csv",index=False)
    else:
        df=run_power(rng,a.reps,a.surr,qs,[float(x) for x in a.sigma.split(',')]); df.to_csv(f"{a.out_dir}/registered_test_power_sensitivity.csv",index=False)
        print("\nInterpretation: power is a sensitivity surface, not an empirical point estimate.")
        print("SSF fixes total smooth variance; q specifies how much of that smooth variance is E2-coupled.")
    print(f"[{time.time()-t0:.1f}s]")

if __name__=='__main__': main()
