#!/usr/bin/env python3
"""Chunked submission-grade H4 theoretical frontier at 500 Monte Carlo replicates."""
import argparse,os,time,warnings
warnings.filterwarnings('ignore')
import numpy as np,pandas as pd
from fastlrt import lrt_random_slope_fast
from calibration_fidelity_aliasing import simulate,within_r,even_days
DENS=[2,3,5,7,10,14,21,28]; RELS=[.55,.60,.65,.70,.75,.80,.85]
def cell(rng,d,rel,sb,n):
    both=c1=c2=valid=0
    for _ in range(n):
        Y,x,_=simulate(rng,even_days(d,2),rel,sb); prop=float(np.nanmean(np.abs(within_r(Y,x))>.20)); p=lrt_random_slope_fast(Y,x)
        if not np.isfinite(p): continue
        valid+=1; a=prop>.50; b=p<.05; c1+=a; c2+=b; both+=(a and b)
    den=max(valid,1); return both/den,c1/den,c2/den,valid
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('sigma_b',type=float); ap.add_argument('chunk',type=int); ap.add_argument('n_chunks',type=int); ap.add_argument('--n-sims',type=int,default=500); ap.add_argument('--out-dir',default='./results/chunks'); ap.add_argument('--seed',type=int,default=7000); a=ap.parse_args(); os.makedirs(a.out_dir,exist_ok=True)
    mine=[(d,r) for d in DENS for r in RELS][a.chunk::a.n_chunks]; rng=np.random.default_rng(a.seed+a.chunk+int(round(a.sigma_b*1000))); rows=[]; t0=time.time()
    for d,r in mine:
        rec,s,l,nvalid=cell(rng,d,r,a.sigma_b,a.n_sims); rows.append(dict(sigma_b=a.sigma_b,obs_per_cycle=d,reliability=r,n_cycles=2,total_obs_per_subj=d*2,recovery=rec,crit1_effectsize=s,crit2_lrt=l,n_sims=a.n_sims,n_valid=nvalid)); print(f'sb={a.sigma_b} d={d} rel={r:.2f} -> rec={rec:.3f}',flush=True)
    out=f'{a.out_dir}/h4_500_sb{a.sigma_b}_c{a.chunk}.csv'; pd.DataFrame(rows).to_csv(out,index=False); print(f'[done] {out} in {time.time()-t0:.0f}s')
if __name__=='__main__': main()
