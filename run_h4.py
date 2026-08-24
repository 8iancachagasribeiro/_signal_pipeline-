#!/usr/bin/env python3
"""H4 theoretical design frontier using the validated fast LRT."""
import warnings,time
warnings.filterwarnings('ignore')
import numpy as np,pandas as pd
import h4_frontier as H
from fastlrt import lrt_random_slope_fast
N_SUBJ=H.N_SUBJ; DENSITIES=[2,3,5,7,10,14,21,28]; RELIABILITIES=[.55,.60,.65,.70,.75,.80,.85]; SIGMAS=[.05,.10,.15,.20]; N_CYCLES=2; N_SIMS=500; SEED=7

def simulate_matrix(rng,n_subj,obs_per_cycle,n_cycles,reliability,sigma_b,balance_offset=0.):
    days=np.concatenate([np.linspace(0,H.CYCLE_LEN,obs_per_cycle,endpoint=False)+c*H.CYCLE_LEN for c in range(n_cycles)]); x=H.e2(days); b=rng.normal(H.DA_OPT-H.K_GAIN*H._E2_MEAN+balance_offset,sigma_b,size=n_subj); signal=H.inverted_u(b[:,None]+H.K_GAIN*x[None,:]); y_true=signal+rng.normal(0,H.SIGMA_STATE,size=signal.shape); sd=np.sqrt(max(y_true.var(),1e-12)*(1-reliability)/reliability); return y_true+rng.normal(0,sd,size=y_true.shape),y_true,x

def within_r_matrix(Y,x):
    xc=x-x.mean(); Yc=Y-Y.mean(axis=1,keepdims=True); den=np.sqrt((Yc**2).sum(axis=1)*(xc**2).sum()); den=np.where(den<1e-12,np.nan,den); return (Yc@xc)/den

def eval_cell(args):
    sb,d,rel,seed,n_sims=args; rng=np.random.default_rng(seed); c1=c2=both=valid=0; grp=[]; props=[]
    for _ in range(n_sims):
        Yo,_,x=simulate_matrix(rng,N_SUBJ,d,N_CYCLES,rel,sb); r_i=within_r_matrix(Yo,x); prop=float(np.nanmean(np.abs(r_i)>.20)); p=lrt_random_slope_fast(Yo,x)
        if not np.isfinite(p): continue
        valid+=1; a=prop>.50; b=p<.05; c1+=a; c2+=b; both+=(a and b); props.append(prop); grp.append(float(np.corrcoef(np.tile(x,N_SUBJ),Yo.ravel())[0,1]))
    den=max(valid,1); return dict(sigma_b=sb,obs_per_cycle=d,n_cycles=N_CYCLES,total_obs_per_subj=d*N_CYCLES,reliability=rel,recovery=both/den,crit1_effectsize=c1/den,crit2_lrt=c2/den,mean_group_r=float(np.mean(grp)) if grp else np.nan,mean_prop_detected=float(np.mean(props)) if props else np.nan,n_valid=valid)

def main():
    import argparse,os
    from multiprocessing import Pool
    ap=argparse.ArgumentParser(); ap.add_argument('sigma_b',nargs='?',type=float,default=None); ap.add_argument('--n-sims',type=int,default=N_SIMS); ap.add_argument('--out-dir',default='./results'); ap.add_argument('--workers',type=int,default=4); a=ap.parse_args(); os.makedirs(a.out_dir,exist_ok=True)
    sigmas=[a.sigma_b] if a.sigma_b is not None else SIGMAS; jobs=[]; k=0
    for sb in sigmas:
        for d in DENSITIES:
            for rel in RELIABILITIES:
                k+=1; jobs.append((sb,d,rel,SEED*1000+k,a.n_sims))
    t0=time.time()
    with Pool(processes=a.workers) as pool: rows=pool.map(eval_cell,jobs)
    df=pd.DataFrame(rows); tag=f'_sb{a.sigma_b}' if a.sigma_b is not None else ''; out=f'{a.out_dir}/h4_part{tag}.csv'; df.to_csv(out,index=False); print(f'[done] {len(df)} cells in {time.time()-t0:.0f}s -> {out}')
if __name__=='__main__': main()
