from __future__ import annotations
import argparse,os,time
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np,pandas as pd
import h4_frontier as H
from brm_common import SEEDS,e2_shift
SPAN=90; ALPHA=.05

def corr_batch(XB,Y,M):
    # XB: B,N,T; Y,M: N,T
    m=M.astype(float); n=np.maximum(m.sum(1),1.0); ym=Y*m; sy=ym.sum(1); sy2=(Y*ym).sum(1)
    sx=(XB*m[None]).sum(2); sx2=((XB*XB)*m[None]).sum(2); sxy=(XB*ym[None]).sum(2)
    num=sxy-sx*sy[None]/n[None]; vx=sx2-sx*sx/n[None]; vy=sy2-sy*sy/n
    return num/np.sqrt(np.maximum(vx,1e-12)*np.maximum(vy[None],1e-12))
def ar_noise(rng,N,T,rho):
    z=rng.normal(size=(N,T));out=np.empty_like(z);out[:,0]=z[:,0]
    sc=np.sqrt(max(1-rho*rho,1e-12))
    for t in range(1,T):out[:,t]=rho*out[:,t-1]+sc*z[:,t]
    return out

def one_test(rng,N,nobs,beta_sd,resid_rho,B,beta_mean=0.0):
    days=np.arange(SPAN,dtype=float);offs=rng.uniform(0,H.CYCLE_LEN,N)
    X=np.array([e2_shift(days+o,0) for o in offs]);X=(X-X.mean(1,keepdims=True))/np.maximum(X.std(1,keepdims=True),1e-9)
    beta=np.full(N,beta_mean) if beta_sd==0 else rng.normal(beta_mean,beta_sd,N)
    Y=beta[:,None]*X+ar_noise(rng,N,SPAN,resid_rho)
    M=np.zeros((N,SPAN),bool)
    if nobs>=SPAN:M[:]=True
    else:
        for i in range(N):M[i,rng.choice(SPAN,nobs,replace=False)]=True
    r_obs=corr_batch(X[None],Y,M)[0];S_obs=float(np.std(r_obs,ddof=0))
    XF=np.fft.rfft(X,axis=1);mag=np.abs(XF);nf=XF.shape[1]
    ph=rng.uniform(0,2*np.pi,size=(B,N,nf));ph[:,:,0]=0
    if SPAN%2==0:ph[:,:,-1]=0
    XS=np.fft.irfft(mag[None]*np.exp(1j*ph),n=SPAN,axis=2)
    sn=corr_batch(XS,Y,M).std(1,ddof=0)
    p=(1+int(np.sum(sn>=S_obs)))/(B+1)
    return p,S_obs,float(np.median(sn))
def run_cell(job):
    kind,params,R,B,seed=job;rng=np.random.default_rng(seed);rej=0;sobs=[];snull=[]
    for _ in range(R):
        p,s,n=one_test(rng,params['N'],params['nobs'],params['beta_sd'],params['rho'],B,params.get('beta_mean',0.0));rej+=p<ALPHA;sobs.append(s);snull.append(n)
    rate=rej/R;mcse=np.sqrt(rate*(1-rate)/R)
    return {**params,'analysis':kind,'rejection_rate':rate,'mcse':mcse,'mean_SD_obs':float(np.mean(sobs)),'mean_null_median':float(np.mean(snull)),'replications':R,'surrogates_per_test':B,'seed':seed}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reps',type=int,default=1000);ap.add_argument('--surrogates',type=int,default=199);ap.add_argument('--workers',type=int,default=max(1,min(8,os.cpu_count() or 2)));ap.add_argument('--out-dir',default='./results')
    a=ap.parse_args();os.makedirs(a.out_dir,exist_ok=True);jobs=[]
    # Calibration mirrors audited mcPHASES scale: 41 participants, 85 observed days in a 90-day calendar.
    for i,rho in enumerate((.2,.5)):
        jobs.append(('calibration',dict(N=41,nobs=85,beta_sd=0.0,beta_mean=0.0,rho=rho),a.reps,a.surrogates,SEEDS['surrogate_calibration']+i*1000))
        # Boundary condition: homogeneous nonzero coupling, still no slope heterogeneity.
        jobs.append(('calibration_common_slope',dict(N=41,nobs=85,beta_sd=0.0,beta_mean=0.2,rho=rho),a.reps,a.surrogates,SEEDS['surrogate_calibration']+1000000+i*1000))
    q=0
    for sd in (.1,.2,.3):
      for N in (20,40,60):
       for nobs in (30,60,90):
        q+=1;jobs.append(('power',dict(N=N,nobs=nobs,beta_sd=sd,beta_mean=0.0,rho=.2),a.reps,a.surrogates,SEEDS['surrogate_power']+q*1009))
    rows=[];t0=time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        fut=[ex.submit(run_cell,j) for j in jobs]
        for k,f in enumerate(as_completed(fut),1):
            rows.append(f.result());print(f'[surrogate] {k}/{len(jobs)} cells',flush=True)
    df=pd.DataFrame(rows);path=os.path.join(a.out_dir,'brm_surrogate_calibration_power.csv');df.to_csv(path,index=False);print(f'[saved] {path}; elapsed={time.time()-t0:.1f}s')
if __name__=='__main__':main()
