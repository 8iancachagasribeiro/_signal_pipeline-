#!/usr/bin/env python3
"""Generative H4 mechanism and theoretical design frontier."""
import argparse,warnings,os
import numpy as np,pandas as pd
from scipy import stats
warnings.filterwarnings('ignore')
CYCLE_LEN=28.; DA_OPT=.50; U_WIDTH=.35; K_GAIN=.15; SIGMA_STATE=.085; N_SUBJ=39
_grid=np.arange(0,CYCLE_LEN,.25)
def _e2_raw(d): return .15+np.exp(-((d-13.)**2)/(2*2.**2))+.55*np.exp(-((d-21.)**2)/(2*3.5**2))
_E2_MAX=_e2_raw(_grid).max(); _E2_MEAN=(_e2_raw(_grid)/_E2_MAX).mean()
def e2(days): return _e2_raw(np.asarray(days,float)%CYCLE_LEN)/_E2_MAX
def inverted_u(da): return np.exp(-((da-DA_OPT)**2)/(2*U_WIDTH**2))
def simulate_study(rng,n_subj,obs_per_cycle,n_cycles,reliability,sigma_b,balance_offset=0.):
    days=np.concatenate([np.linspace(0,CYCLE_LEN,obs_per_cycle,endpoint=False)+c*CYCLE_LEN for c in range(n_cycles)]); x=e2(days); n_obs=len(days)
    b=rng.normal(DA_OPT-K_GAIN*_E2_MEAN+balance_offset,sigma_b,size=n_subj); signal=inverted_u(b[:,None]+K_GAIN*x[None,:]); y_true=signal+rng.normal(0,SIGMA_STATE,size=signal.shape)
    sd_err=np.sqrt(max(y_true.var(),1e-12)*(1-reliability)/reliability); y_obs=y_true+rng.normal(0,sd_err,size=y_true.shape)
    return pd.DataFrame({'subj':np.repeat(np.arange(n_subj),n_obs),'e2':np.tile(x,n_subj),'y_true':y_true.ravel(),'y_obs':y_obs.ravel()})
def within_person_r(df,col='y_obs'):
    out=[]
    for _,g in df.groupby('subj',sort=False): out.append(0. if g['e2'].std()<1e-9 or g[col].std()<1e-9 else np.corrcoef(g['e2'],g[col])[0,1])
    return np.asarray(out)
def group_effect(df,col='y_obs'): return 0. if df['e2'].std()<1e-9 or df[col].std()<1e-9 else float(np.corrcoef(df['e2'],df[col])[0,1])
def lrt_random_slope(df,col='y_obs'):
    import statsmodels.formula.api as smf
    try:
        m0=smf.mixedlm(f'{col} ~ e2',df,groups=df['subj']).fit(reml=False); m1=smf.mixedlm(f'{col} ~ e2',df,groups=df['subj'],re_formula='~e2').fit(reml=False); stat=2*(m1.llf-m0.llf)
        if not np.isfinite(stat) or stat<0: return np.nan
        return float(.5*stats.chi2.sf(stat,1)+.5*stats.chi2.sf(stat,2))
    except Exception: return np.nan
def evaluate_cell(rng,n_sims,obs_per_cycle,n_cycles,reliability,sigma_b,balance_offset=0.,n_subj=N_SUBJ):
    c1=c2=both=ok=0; grp=[]; propr=[]
    for _ in range(n_sims):
        df=simulate_study(rng,n_subj,obs_per_cycle,n_cycles,reliability,sigma_b,balance_offset); r_i=within_person_r(df); prop=float(np.mean(np.abs(r_i)>.20)); p=lrt_random_slope(df)
        if np.isnan(p): continue
        ok+=1; a=prop>.50; b=p<.05; c1+=a; c2+=b; both+=(a and b); grp.append(group_effect(df)); propr.append(prop)
    if ok==0: return dict(recovery=np.nan,crit1=np.nan,crit2=np.nan,mean_group_r=np.nan,mean_prop=np.nan,n_ok=0)
    return dict(recovery=both/ok,crit1=c1/ok,crit2=c2/ok,mean_group_r=float(np.mean(grp)),mean_prop=float(np.mean(propr)),n_ok=ok)
def validate_masking(seed=0):
    rng=np.random.default_rng(seed); print('='*74); print('STEP A - masking sanity check'); print('='*74)
    print(f"{'sigma_b':>8} {'group r':>9} {'|group g|':>10} {'med |r_i|':>10} {'% pos':>8}")
    for sb in [.05,.10,.15,.20]:
        gr=[]; mr=[]; ps=[]
        for _ in range(200):
            df=simulate_study(rng,N_SUBJ,28,2,.73,sb); r=within_person_r(df,'y_true'); gr.append(group_effect(df,'y_true')); mr.append(np.median(np.abs(r))); ps.append(np.mean(r>0))
        rr=float(np.mean(gr)); gg=abs(2*rr/np.sqrt(max(1-rr**2,1e-9))); print(f'{sb:>8.2f} {rr:>9.3f} {gg:>10.3f} {np.mean(mr):>10.3f} {100*np.mean(ps):>7.1f}%')
    print('\nSTEP B - directional imbalance boundary')
    for off in [0.,.05,.10,.15,.20,.30]:
        gr=[]; ps=[]
        for _ in range(200):
            df=simulate_study(rng,N_SUBJ,28,2,.73,.10,balance_offset=off); r=within_person_r(df,'y_true'); gr.append(group_effect(df,'y_true')); ps.append(np.mean(r>0))
        rr=float(np.mean(gr)); gg=abs(2*rr/np.sqrt(max(1-rr**2,1e-9))); print(f'offset={off:.2f} |g|={gg:.3f} pos={100*np.mean(ps):.1f}%')
def run_frontier(seed,n_sims,densities,reliabilities,sigmas,n_cycles):
    rng=np.random.default_rng(seed); rows=[]
    for sb in sigmas:
        for d in densities:
            for r in reliabilities:
                res=evaluate_cell(rng,n_sims,d,n_cycles,r,sb); rows.append(dict(sigma_b=sb,obs_per_cycle=d,n_cycles=n_cycles,total_obs=d*n_cycles,reliability=r,**res)); print(f'sb={sb:.2f} d={d} rel={r:.2f} recovery={res["recovery"]:.2f}',flush=True)
    return pd.DataFrame(rows)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--validate-only',action='store_true'); ap.add_argument('--quick',action='store_true'); ap.add_argument('--n-sims',type=int,default=100); ap.add_argument('--seed',type=int,default=7); ap.add_argument('--n-cycles',type=int,default=2); ap.add_argument('--out',default='./results/h4_frontier_results.csv'); a=ap.parse_args(); validate_masking(a.seed)
    if a.validate_only: return
    if a.quick: densities=[2,5,14]; reliabilities=[.55,.73,.85]; sigmas=[.10]; n_sims=20
    else: densities=[2,3,5,7,10,14,21,28]; reliabilities=[.55,.60,.65,.70,.75,.80,.85]; sigmas=[.05,.10,.15,.20]; n_sims=a.n_sims
    df=run_frontier(a.seed,n_sims,densities,reliabilities,sigmas,a.n_cycles); os.makedirs(os.path.dirname(a.out) or '.',exist_ok=True); df.to_csv(a.out,index=False); print(f'[saved] {a.out}')
if __name__=='__main__': main()
