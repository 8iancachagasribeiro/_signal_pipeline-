#!/usr/bin/env python3
"""Analytical attenuation sensitivity and design-budget sensitivity."""
from __future__ import annotations
import argparse, os
import numpy as np, pandas as pd
SSF_PREDICTOR=.469; SSF_OUTCOME=.323

def table_9():
    rows=[]; print(f"{'f':>6} {'rel predictor':>14} {'rel outcome':>12} {'attenuation':>12} {'r_true for .20':>15}")
    for f in (0.,.25,.50,.75,1.):
        rx=SSF_PREDICTOR+f*(1-SSF_PREDICTOR); ry=SSF_OUTCOME+f*(1-SSF_OUTCOME); att=np.sqrt(rx*ry)
        rows.append(dict(f=f,reliability_predictor=rx,reliability_outcome=ry,attenuation=att,r_true_for_r_obs_020=.20/att))
        print(f"{f:>6.2f} {rx:>14.3f} {ry:>12.3f} {att:>12.3f} {.20/att:>15.3f}")
    print('\nThis table is an assumption analysis, not an estimate of psychometric reliability.')
    return pd.DataFrame(rows)

def table_17(n_sims,q_grid,seed):
    import h4_v2 as V
    return V.budget_sensitivity(np.random.default_rng(seed),n_sims,q_grid)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--table',type=int,choices=[9,17],default=9); ap.add_argument('--n-sims',type=int,default=200); ap.add_argument('--q',default='0.10,0.25,0.50,0.75,1.00'); ap.add_argument('--seed',type=int,default=2026); ap.add_argument('--out-dir',default='./results')
    a=ap.parse_args(); os.makedirs(a.out_dir,exist_ok=True)
    if a.table==9: df=table_9(); path=f'{a.out_dir}/table09_attenuation_sensitivity.csv'
    else:
        df=table_17(a.n_sims,[float(x) for x in a.q.split(',')],a.seed); path=f'{a.out_dir}/table17_budget_sensitivity.csv'
        print("\nNo single row is 'the empirical power': q is not identified by SSF alone.")
    df.to_csv(path,index=False); print(f'[saved] {path}')
if __name__=='__main__': main()
