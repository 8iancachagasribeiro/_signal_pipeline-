#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np, pandas as pd

def near(x,y,tol=5e-4):
    if not np.isfinite(x) or abs(x-y)>tol: raise AssertionError(f'{x} != {y} within {tol}')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',default='.');a=ap.parse_args();root=Path(a.root).resolve();res=root/'results'
    rob=pd.read_csv(res/'brm_robustness_metrics.csv'); samp=pd.read_csv(res/'brm_sampling_design.csv'); ssf=pd.read_csv(res/'brm_ssf_benchmark.csv'); sur=pd.read_csv(res/'brm_surrogate_calibration_power.csv')
    assert len(rob)==400 and (rob.replications==1000).all()
    assert len(samp)==384 and (samp.replications==1000).all()
    assert len(ssf)==576 and (ssf.replications==1000).all()
    assert len(sur)==31 and (sur.replications==1000).all() and (sur.surrogates_per_test==199).all()
    ref=rob[(rob.factor=='reference')&(rob.level=='reference')&(rob.n_cycles==2)].set_index('n_obs')
    near(ref.loc[20,'fidelity'],.698806,1e-6); near(ref.loc[56,'fidelity'],.812056,1e-6)
    sig=rob[(rob.factor=='functional_form')&(rob.level=='sigmoid')&(rob.n_cycles==2)].set_index('n_obs')
    near(sig.loc[56,'fidelity'],.603824,1e-6)
    pt=samp[(samp.schedule=='phase_targeted')&(samp.obs_per_cycle==10)&(samp.missingness=='none')].sort_values('peak_sigma_days')
    exp=[.706202,.704328,.700096,.705570]
    for x,y in zip(pt.fidelity,exp): near(x,y,1e-6)
    agg=ssf.groupby('method').mse.mean(); near(agg['ssf_020'],.020135,1e-6); near(agg['acf_linear'],.021479,1e-6)
    cal=sur[sur.analysis=='calibration'].set_index('rho'); near(cal.loc[.2,'rejection_rate'],.057,1e-12); near(cal.loc[.5,'rejection_rate'],.036,1e-12)
    p=sur[(sur.analysis=='power')&(sur.beta_sd==.2)&(sur.N==20)].set_index('nobs');
    near(p.loc[30,'rejection_rate'],.597,1e-12); near(p.loc[60,'rejection_rate'],.863,1e-12); near(p.loc[90,'rejection_rate'],.948,1e-12)
    print('VALIDATION PASS: canonical row counts, replication counts, and manuscript anchor values match.')

if __name__=='__main__': main()
