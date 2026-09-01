from __future__ import annotations
import json, math, os
import numpy as np
import h4_frontier as H

REL_REF=.73
BASE_B=H.DA_OPT-H.K_GAIN*H._E2_MEAN
DENSE_X=H.e2(np.linspace(0,H.CYCLE_LEN,400,endpoint=False))

SEEDS={
    'robustness': 26083101,
    'ssf_benchmark': 26083102,
    'sampling_design': 26083103,
    'surrogate_calibration': 26083104,
    'surrogate_power': 26083105,
}

def response(da, form='inverted_u'):
    da=np.asarray(da,float)
    if form=='inverted_u': return H.inverted_u(da)
    if form=='sigmoid': return 1/(1+np.exp(-10*(da-H.DA_OPT)))
    if form=='asym_linear':
        z=da-H.DA_OPT
        return np.where(z<0,1.10*z,-0.65*z)
    if form=='piecewise':
        z=da-H.DA_OPT
        return np.where(z<-0.08,0.70*z,np.where(z<0.10,1.0*z,-0.50*z+0.15))
    raise ValueError(form)

def e2_shift(days, shift=0.0):
    d=np.asarray(days,float)%H.CYCLE_LEN
    raw=.15+np.exp(-((d-(13.+shift))**2)/(2*2.**2))+.55*np.exp(-((d-21.)**2)/(2*3.5**2))
    return raw/H._E2_MAX

def corr_rows(Y,x):
    Y=np.asarray(Y,float); x=np.asarray(x,float)
    if x.ndim==1: x=np.broadcast_to(x,Y.shape)
    m=np.isfinite(Y)&np.isfinite(x)
    n=m.sum(-1).astype(float)
    n_safe=np.maximum(n,1)
    yy=np.where(m,Y,0.0); xx=np.where(m,x,0.0)
    sy=yy.sum(-1); sx=xx.sum(-1); sy2=(yy*yy).sum(-1); sx2=(xx*xx).sum(-1); sxy=(yy*xx).sum(-1)
    num=sxy-sx*sy/n_safe
    vy=sy2-sy*sy/n_safe; vx=sx2-sx*sx/n_safe
    r=num/np.sqrt(np.maximum(vy,1e-14)*np.maximum(vx,1e-14))
    r=np.where((n>=4)&(vy>1e-12)&(vx>1e-12),r,np.nan)
    return np.clip(r,-1,1)

def corr_pairwise_rows(A,B):
    return corr_rows(A,B)

def _build_truth_lookup(form):
    bg=np.linspace(.05,.95,6001)
    sig=response(bg[:,None]+H.K_GAIN*DENSE_X[None,:],form)
    rr=corr_rows(sig,np.broadcast_to(DENSE_X,sig.shape))
    return bg,rr
_TRUTH={f:_build_truth_lookup(f) for f in ('inverted_u','sigmoid','asym_linear','piecewise')}

def true_coupling(b,K=None,form='inverted_u'):
    b=np.asarray(b,float)
    if K is None or (np.ndim(K)==0 and abs(float(K)-H.K_GAIN)<1e-12):
        bg,rr=_TRUTH[form]
        return np.interp(b,bg,rr,left=rr[0],right=rr[-1])
    K=np.asarray(K,float)
    flatb=b.ravel(); flatk=np.broadcast_to(K,b.shape).ravel(); out=np.empty_like(flatb)
    chunk=4000
    for a in range(0,len(flatb),chunk):
        bb=flatb[a:a+chunk,None]; kk=flatk[a:a+chunk,None]
        sig=response(bb+kk*DENSE_X[None,:],form)
        out[a:a+chunk]=corr_rows(sig,np.broadcast_to(DENSE_X,sig.shape))
    return out.reshape(b.shape)

def ar1_noise(rng,shape,rho,sd):
    z=rng.normal(size=shape)
    if rho==0: return z*sd
    out=np.empty(shape,float); out[...,0]=z[...,0]*sd
    innov=sd*np.sqrt(max(1-rho*rho,1e-12))
    for t in range(1,shape[-1]): out[...,t]=rho*out[...,t-1]+innov*z[...,t]
    return out

def rep_metrics(est,true):
    est=np.asarray(est,float); true=np.asarray(true,float)
    R=est.shape[0]
    fidelity=np.full(R,np.nan); rmse=np.full(R,np.nan); bias=np.full(R,np.nan); sign=np.full(R,np.nan); icc=np.full(R,np.nan)
    for r in range(R):
        m=np.isfinite(est[r])&np.isfinite(true[r])
        a=est[r,m]; b=true[r,m]
        if len(a)<4: continue
        if np.std(a)>1e-12 and np.std(b)>1e-12: fidelity[r]=np.corrcoef(a,b)[0,1]
        d=a-b; rmse[r]=np.sqrt(np.mean(d*d)); bias[r]=np.mean(d); sign[r]=np.mean(np.sign(a)==np.sign(b))
        X=np.column_stack([b,a]); n,k=X.shape
        gm=X.mean(); subj=X.mean(1); rat=X.mean(0)
        ssr=k*np.sum((subj-gm)**2); sse=np.sum((X-subj[:,None]-rat[None,:]+gm)**2)
        msr=ssr/max(n-1,1); mse=sse/max((n-1)*(k-1),1)
        den=msr+(k-1)*mse
        icc[r]=(msr-mse)/den if den>1e-12 else np.nan
    return {'fidelity':fidelity,'rmse':rmse,'bias':bias,'icc':icc,'sign_accuracy':sign}

def summarize_metrics(metrics,R):
    out={}
    for name,v in metrics.items():
        v=np.asarray(v,float); f=v[np.isfinite(v)]
        out[name]=float(np.mean(f)) if len(f) else np.nan
        out[name+'_sd']=float(np.std(f,ddof=1)) if len(f)>1 else np.nan
        out[name+'_mcse']=float(np.std(f,ddof=1)/np.sqrt(len(f))) if len(f)>1 else np.nan
        out[name+'_nvalid']=int(len(f))
    out['replications']=int(R)
    return out

def write_seed_map(path):
    with open(path,'w',encoding='utf-8') as f: json.dump(SEEDS,f,indent=2,sort_keys=True)
