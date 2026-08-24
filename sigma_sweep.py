#!/usr/bin/env python3
"""Sensitivity analysis over heterogeneity and the unidentified E2-coupled smooth share q."""
import argparse,os
import numpy as np,pandas as pd
from registered_test_power import simulate_study,surrogate_test,SSF_X,SSF_Y
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--q',default='0.05,0.10,0.25,0.50,0.75,1.00'); ap.add_argument('--sigma',default='0.04,0.075,0.12,0.20'); ap.add_argument('--reps',type=int,default=100); ap.add_argument('--surr',type=int,default=199); ap.add_argument('--seed',type=int,default=41); ap.add_argument('--out-dir',default='./results'); a=ap.parse_args(); os.makedirs(a.out_dir,exist_ok=True)
    rng=np.random.default_rng(a.seed); rows=[]
    for q in [float(v) for v in a.q.split(',')]:
        for sb in [float(v) for v in a.sigma.split(',')]:
            rej=0; spread=[]
            for _ in range(a.reps):
                rec=simulate_study(rng,sb,coupled_fraction=q); p,s,_=surrogate_test(rec,rng,B=a.surr); rej+=p<.05; spread.append(s)
            rows.append(dict(q=q,sigma_b=sb,power=rej/a.reps,mean_SD_ri=np.mean(spread),target_ssf_predictor=SSF_X,target_ssf_outcome=SSF_Y,reps=a.reps,surrogates=a.surr)); print(f'q={q:.2f} sigma_b={sb:.3f} -> power={rej/a.reps:.3f}')
    path=f'{a.out_dir}/q_sigma_sensitivity.csv'; pd.DataFrame(rows).to_csv(path,index=False); print(f'[saved] {path}')
if __name__=='__main__': main()
