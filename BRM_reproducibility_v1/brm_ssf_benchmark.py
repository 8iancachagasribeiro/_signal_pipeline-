from __future__ import annotations
import argparse, os, time, math
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np, pandas as pd
from scipy.interpolate import UnivariateSpline
from statsmodels.tsa.statespace.structural import UnobservedComponents
from brm_common import SEEDS
LN2=np.log(2.0)
LENGTHS=(30,60,90,120,160,200); TRUE_FRACS=(.2,.4,.6,.8); SHAPES=('sinusoid','two_peak','pulse')

def signal_shape(name,n):
    t=np.arange(n,dtype=float); d=t%28
    if name=='sinusoid': s=np.sin(2*np.pi*t/28)
    elif name=='two_peak': s=np.exp(-((d-13.)**2)/8)+.55*np.exp(-((d-21.)**2)/24)
    elif name=='pulse': s=(np.minimum(d,28-d)<2.5).astype(float)
    else: raise ValueError(name)
    s=s-s.mean(); sd=s.std(); return s/(sd if sd>1e-12 else 1)

def ssf_spectral_raw(y,cut=.25):
    y=np.asarray(y,float); n=len(y)
    if n<8 or np.std(y)<1e-12: return np.nan
    z=y-y.mean(); P=np.abs(np.fft.rfft(z))**2/n; f=np.fft.rfftfreq(n,d=1.0); P=P[1:]; f=f[1:]
    hi=f>cut
    if hi.sum()<4 or P.sum()<=0: return np.nan
    noise=np.median(P[hi])/LN2
    return 1-noise*len(P)/P.sum()

def adaptive_cut(y):
    y=np.asarray(y,float); n=len(y); z=y-y.mean(); P=np.abs(np.fft.rfft(z))**2/n; f=np.fft.rfftfreq(n,d=1.0); P=P[1:];f=f[1:]
    candidates=[.15,.20,.25,.30,.35]
    best=(1e9,.25)
    for c in candidates:
        q=P[f>c]
        if len(q)<8: continue
        h=len(q)//2; a=np.median(q[:h])+1e-12; b=np.median(q[h:])+1e-12
        score=abs(np.log(a/b))+.08*c
        if score<best[0]: best=(score,c)
    return best[1]

def acf_linear(y):
    y=np.asarray(y,float); n=len(y)
    if n<8 or np.std(y)<1e-12:return np.nan
    z=y-y.mean(); v=np.dot(z,z)/n
    r1=np.dot(z[:-1],z[1:])/n/v; r2=np.dot(z[:-2],z[2:])/n/v
    return 2*r1-r2

def state_space_frac(y):
    y=np.asarray(y,float); n=len(y)
    if n<8 or np.std(y)<1e-12:return np.nan
    dy=np.diff(y); vd=np.var(dy)
    if len(dy)>2:
        cov1=np.mean((dy[:-1]-dy[:-1].mean())*(dy[1:]-dy[1:].mean()))
    else: cov1=0
    obs=max(-cov1,1e-6*np.var(y)); level=max(vd-2*obs,1e-6*np.var(y))
    try:
        m=UnobservedComponents(y,level=True,irregular=True,stochastic_level=True)
        r=m.smooth([obs,level]); s=np.asarray(r.smoothed_state[0])
        return float(np.var(s)/np.var(y))
    except Exception: return np.nan

def _haar_denoise(y):
    y=np.asarray(y,float); n0=len(y); m=1
    while m<n0: m*=2
    if m>n0: x=np.pad(y,(0,m-n0),mode='reflect')
    else: x=y.copy()
    coeff=[]; a=x.copy()
    while len(a)>=4:
        av=(a[0::2]+a[1::2])/np.sqrt(2); d=(a[0::2]-a[1::2])/np.sqrt(2); coeff.append(d); a=av
    finest=coeff[0] if coeff else np.array([0.])
    sigma=np.median(np.abs(finest-np.median(finest)))/.67448975 if len(finest)>2 else np.std(finest)
    thr=max(sigma,1e-12)*np.sqrt(2*np.log(max(m,2)))
    coeff=[np.sign(d)*np.maximum(np.abs(d)-thr,0) for d in coeff]
    for d in coeff[::-1]:
        out=np.empty(len(d)*2); out[0::2]=(a+d)/np.sqrt(2); out[1::2]=(a-d)/np.sqrt(2); a=out
    return a[:n0]
def wavelet_frac(y):
    try: s=_haar_denoise(y); return float(np.var(s)/np.var(y)) if np.var(y)>1e-12 else np.nan
    except Exception:return np.nan

def spline_frac(y):
    """Penalized cubic spline with smoothing selected by blocked 3-fold CV.

    scipy.interpolate.UnivariateSpline is used for speed. Candidate smoothing budgets
    are scaled by n*Var(y), making the grid invariant to the measurement scale.
    """
    y=np.asarray(y,float); n=len(y)
    if n<12 or np.std(y)<1e-12:return np.nan
    x=np.arange(n,dtype=float); vy=float(np.var(y))
    scales=(.10,.30,1.0,3.0,10.0)
    # Interleaved folds preserve coverage over the full time span while withholding points.
    folds=np.arange(n)%3; best=None
    try:
        for sc in scales:
            errs=[]
            for f in range(3):
                tr=folds!=f; te=~tr
                spl=UnivariateSpline(x[tr],y[tr],k=3,s=float(sc*tr.sum()*vy))
                pr=spl(x[te]); errs.append(float(np.mean((y[te]-pr)**2)))
            score=float(np.mean(errs))
            if best is None or score<best[0]: best=(score,sc)
        sc=best[1]; spl=UnivariateSpline(x,y,k=3,s=float(sc*n*vy)); sm=spl(x)
        return float(np.var(sm)/vy)
    except Exception:return np.nan

def _cell(args):
    shape,truef,n,R,seed=args; rng=np.random.default_rng(seed); s=signal_shape(shape,n); noise_sd=np.sqrt((1-truef)/truef)
    methods=['ssf_020','ssf_025','ssf_030','ssf_adaptive','state_space','haar_wavelet','spline_gcv','acf_linear']
    vals={m:[] for m in methods}; fail={m:0 for m in methods}
    chosen=[]
    for _ in range(R):
        y=s+rng.normal(0,noise_sd,n); cut=adaptive_cut(y); chosen.append(cut)
        est={
            'ssf_020':ssf_spectral_raw(y,.20),'ssf_025':ssf_spectral_raw(y,.25),'ssf_030':ssf_spectral_raw(y,.30),
            'ssf_adaptive':ssf_spectral_raw(y,cut),'state_space':state_space_frac(y),'haar_wavelet':wavelet_frac(y),
            'spline_gcv':spline_frac(y),'acf_linear':acf_linear(y)}
        for m,v in est.items():
            bad=(not np.isfinite(v)) or (v<0) or (v>1)
            if bad: fail[m]+=1
            else: vals[m].append(v)
    rows=[]
    for m in methods:
        a=np.asarray(vals[m],float); err=a-truef
        rows.append(dict(shape=shape,true_fraction=truef,n=n,method=m,bias=float(np.mean(err)) if len(a) else np.nan,
                         mse=float(np.mean(err*err)) if len(a) else np.nan,stability_var=float(np.var(a,ddof=1)) if len(a)>1 else np.nan,
                         failure_rate=fail[m]/R,mean_estimate=float(np.mean(a)) if len(a) else np.nan,n_valid=len(a),replications=R,
                         bias_mcse=float(np.std(err,ddof=1)/np.sqrt(len(a))) if len(a)>1 else np.nan,seed=seed,
                         adaptive_cut_median=float(np.median(chosen)) if m=='ssf_adaptive' else np.nan))
    return rows

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--reps',type=int,default=1000); ap.add_argument('--workers',type=int,default=max(1,min(8,os.cpu_count() or 2))); ap.add_argument('--out-dir',default='./results'); ap.add_argument('--seed',type=int,default=SEEDS['ssf_benchmark']); ap.add_argument('--shape-only',default=None); ap.add_argument('--tag',default='')
    a=ap.parse_args(); os.makedirs(a.out_dir,exist_ok=True); jobs=[]; idx=0
    shapes=[a.shape_only] if a.shape_only else list(SHAPES)
    for sh in shapes:
      for tf in TRUE_FRACS:
       for n in LENGTHS:
        idx+=1; jobs.append((sh,tf,n,a.reps,a.seed+idx*997))
    rows=[]; t0=time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
      fut=[ex.submit(_cell,j) for j in jobs]
      for k,f in enumerate(as_completed(fut),1):
        rows.extend(f.result())
        if k%8==0: print(f'[ssf] {k}/{len(jobs)} cells',flush=True)
    df=pd.DataFrame(rows); path=os.path.join(a.out_dir,f'brm_ssf_benchmark{a.tag}.csv'); df.to_csv(path,index=False)
    print(f'[saved] {path}; rows={len(df)}; elapsed={time.time()-t0:.1f}s')
if __name__=='__main__':main()
