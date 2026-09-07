from __future__ import annotations
import argparse, os, time
import numpy as np,pandas as pd
import h4_frontier as H
from brm_common import *
KS=tuple(range(3,11)); PEAK_SIGMAS=(0.,1.,2.,3.); MISS=('none','mcar10','mcar25','weekend'); SCHEDULES=('uniform','phase_targeted','adaptive'); N_CYCLES=2
PHASE_TARGETED={3:[3,13,21],4:[3,9,13,21],5:[3,9,13,17,21],6:[2,7,11,13,18,23],7:[2,6,9,13,17,21,25],8:[1,5,9,12,14,18,22,26],9:[1,4,7,10,13,16,19,22,25],10:list(np.rint(np.linspace(1,27,10)).astype(int))}

def _unique_nearest(vals,k):
    used=[]
    for v in vals:
        order=np.argsort(np.abs(np.arange(28)-v))
        for d in order:
            if int(d) not in used: used.append(int(d)); break
    if len(used)<k:
        for d in range(28):
            if d not in used: used.append(d)
            if len(used)==k:break
    return np.array(sorted(used[:k]),float)
def schedule(k,kind):
    if kind=='uniform': return _unique_nearest(np.linspace(0,27,k),k)
    if kind=='phase_targeted': return np.array(PHASE_TARGETED[k],float)
    grid=np.linspace(0,27,2801); vals=e2_shift(grid,0); der=np.abs(np.gradient(vals,grid)); w=der+.02*der.max(); c=np.cumsum(w); c/=c[-1]
    qs=(np.arange(k)+.5)/k; pts=np.interp(qs,c,grid); return _unique_nearest(pts,k)
def all_days(k,kind):
    base=schedule(k,kind); return np.concatenate([base+c*28 for c in range(N_CYCLES)])

def simulate_cell(seed,k,kind,peak_sigma,missing,R=1000,N=H.N_SUBJ):
    rng=np.random.default_rng(seed); days=all_days(k,kind); T=len(days)
    shifts=np.zeros((R,N)) if peak_sigma==0 else rng.normal(0,peak_sigma,size=(R,N))
    dmod=days[None,None,:]%28
    raw=.15+np.exp(-((dmod-(13.+shifts[...,None]))**2)/(2*2.**2))+.55*np.exp(-((dmod-21.)**2)/(2*3.5**2))
    X=raw/H._E2_MAX
    b=rng.normal(BASE_B,.10,size=(R,N)); sig=response(b[...,None]+H.K_GAIN*X,'inverted_u')
    y_true=sig+rng.normal(0,H.SIGMA_STATE,size=sig.shape)
    v=np.var(y_true,axis=(1,2),keepdims=True); sd=np.sqrt(np.maximum(v,1e-12)*(1-REL_REF)/REL_REF); y=y_true+rng.normal(size=y_true.shape)*sd
    if missing=='mcar10': mask=rng.random(y.shape)<.10
    elif missing=='mcar25': mask=rng.random(y.shape)<.25
    elif missing=='weekend':
        wk=(np.rint(days).astype(int)%7>=5); mask=np.broadcast_to(wk,y.shape)
    else: mask=np.zeros(y.shape,bool)
    y=np.where(mask,np.nan,y); Xo=np.where(mask,np.nan,X)
    est=corr_rows(y,Xo); truth=true_coupling(b)
    mets=rep_metrics(est,truth); out=summarize_metrics(mets,R)
    out.update(obs_per_cycle=k,planned_total_obs=T,schedule=kind,peak_sigma_days=peak_sigma,missingness=missing,
               mean_observed=float(np.mean(np.sum(~mask,axis=-1))),median_e2_range=float(np.nanmedian(np.nanmax(Xo,axis=-1)-np.nanmin(Xo,axis=-1))))
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--reps',type=int,default=1000);ap.add_argument('--out-dir',default='./results');ap.add_argument('--seed',type=int,default=SEEDS['sampling_design']);ap.add_argument('--schedule-only',default=None);ap.add_argument('--tag',default='')
    a=ap.parse_args();os.makedirs(a.out_dir,exist_ok=True);rows=[];selected=[a.schedule_only] if a.schedule_only else list(SCHEDULES);total=len(KS)*len(PEAK_SIGMAS)*len(MISS)*len(selected);q=0;t0=time.time()
    for si,kind in enumerate(selected):
      for ps in PEAK_SIGMAS:
       for mi,miss in enumerate(MISS):
        for k in KS:
         seed=a.seed+si*100000+int(ps*10000)+mi*1000+k;rec=simulate_cell(seed,k,kind,ps,miss,a.reps);rec['seed']=seed;rows.append(rec);q+=1
         if q%24==0:print(f'[sampling] {q}/{total} cells',flush=True)
    df=pd.DataFrame(rows);path=os.path.join(a.out_dir,f'brm_sampling_design{a.tag}.csv');df.to_csv(path,index=False);print(f'[saved] {path}; cells={len(df)}; elapsed={time.time()-t0:.1f}s')
if __name__=='__main__':main()
