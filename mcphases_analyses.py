#!/usr/bin/env python3
"""Reproducible empirical analyses for credentialed mcPHASES data.

Key audit corrections:
1) the self-report ordinal map has six levels;
2) SSF is computed on complete regular-grid series, not E3G-paired subsets;
3) phase-randomized surrogates preserve the real calendar rather than compacting gaps;
4) pooled group associations are person-centered;
5) SSF is not treated as classical reliability. Attenuation is reported as a
   sensitivity analysis over the unknown genuine high-frequency fraction f.
"""
from __future__ import annotations
import argparse, os, warnings
import numpy as np, pandas as pd
from ssf_estimators import ssf_spectral, ssf_ar1, ssf_acf_linear, regular_grid
warnings.filterwarnings('ignore')
ORDINAL={'Not at all':0,'Very Low/Little':1,'Low':2,'Moderate':3,'High':4,'Very High':5}
CLASSIFICATION={'fatigue':'BALANCED (confirmatory)','moodswing':'BALANCED (confirmatory)','cramps':'DIRECTIONAL','bloating':'DIRECTIONAL','sorebreasts':'DIRECTIONAL','sleepissue':'AMBIGUOUS (exploratory)','stress':'AMBIGUOUS (exploratory)','appetite':'AMBIGUOUS (exploratory)','foodcravings':'AMBIGUOUS (exploratory)'}
MIN_PAIRED=25

def load(data_dir): return pd.read_csv(f'{data_dir}/hormones_and_selfreport.csv').sort_values(['id','day_in_study'])
def corr(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if np.std(a)<1e-12 or np.std(b)<1e-12: return 0.
    return float(np.corrcoef(a,b)[0,1])
def effect_g(r):
    r=float(np.clip(r,-.99,.99)); return abs(2*r/np.sqrt(1-r*r))
def phase_randomize(x,rng):
    x=np.asarray(x,float); n=len(x); X=np.fft.rfft(x); mag=np.abs(X); ph=rng.uniform(0,2*np.pi,len(X)); ph[0]=0.
    if n%2==0: ph[-1]=0.
    return np.fft.irfft(mag*np.exp(1j*ph),n)
def phase_randomize_indexed(x,index,rng):
    """Phase-randomize on the integer-day calendar, then sample original observed days."""
    x=np.asarray(x,float); idx=np.asarray(index,int)
    if len(x)!=len(idx): raise ValueError('x and index must have equal length')
    order=np.argsort(idx,kind='mergesort'); x,idx=x[order],idx[order]; keep=np.concatenate(([True],np.diff(idx)!=0)); x,idx=x[keep],idx[keep]
    lo,hi=int(idx.min()),int(idx.max()); grid=np.full(hi-lo+1,np.nan); pos=idx-lo; grid[pos]=x; finite=np.isfinite(grid)
    if finite.sum()<4: return x.copy()
    if not finite.all():
        z=np.arange(len(grid)); grid[~finite]=np.interp(z[~finite],z[finite],grid[finite])
    return phase_randomize(grid,rng)[pos]
def paired_series_with_days(df,item,min_n=MIN_PAIRED):
    d=df.dropna(subset=['estrogen',item,'day_in_study']).copy(); d['y']=d[item].map(ORDINAL); d=d.dropna(subset=['y']); out=[]
    for _,g in d.groupby('id'):
        g=g.sort_values('day_in_study')
        if len(g)>=min_n: out.append((g.estrogen.to_numpy(float),g.y.to_numpy(float),g.day_in_study.to_numpy(int)))
    return out
def paired_series(df,item,min_n=MIN_PAIRED): return [(x,y) for x,y,_ in paired_series_with_days(df,item,min_n)]
def pooled_within_corr(per):
    xs=[]; ys=[]
    for rec in per:
        x,y=rec[0],rec[1]; xs.append(x-np.mean(x)); ys.append(y-np.mean(y))
    return corr(np.concatenate(xs),np.concatenate(ys)) if xs else np.nan

def describe(df):
    pairs=paired_series(df,'fatigue'); cnt=pd.Series([len(x) for x,_ in pairs],dtype=float)
    print(f'participants={df.id.nunique()} rows={len(df)} qualifying_fatigue={len(pairs)} median_paired={cnt.median():.0f} IQR={cnt.quantile(.25):.0f}-{cnt.quantile(.75):.0f}')
    print('self-report levels:',sorted(df.fatigue.dropna().unique()))
    print("Preregistered 'energy' confirmatory outcome is absent from mcPHASES and cannot be tested.")

def instrument_ssf(df,data_dir):
    rows=[]
    def report(label,series):
        vals={'ar1':[],'lin':[],'spec':[]}
        for y in series:
            for k,f in [('ar1',ssf_ar1),('lin',ssf_acf_linear),('spec',ssf_spectral)]:
                v=f(y)
                if np.isfinite(v): vals[k].append(v)
        a,l,s=[float(np.median(vals[k])) if vals[k] else np.nan for k in ('ar1','lin','spec')]; rows.append(dict(measure=label,ar1=a,acf_linear=l,spectral=s)); print(f'{label}: spectral SSF={s:.3f}')
    report('estrone-3-glucuronide (PREDICTOR)',[regular_grid(g.estrogen.values,g.day_in_study.values) for _,g in df.groupby('id')])
    for item in ('fatigue','moodswing','cramps','bloating','stress'):
        report(f'{item} [{CLASSIFICATION.get(item,"")}]',[regular_grid(g[item].map(ORDINAL).values,g.day_in_study.values) for _,g in df.groupby('id')])
    try:
        rhr=pd.read_csv(f'{data_dir}/resting_heart_rate.csv').sort_values(['id','day_in_study']); report('resting heart rate (OBJECTIVE)',[regular_grid(g.value.values,g.day_in_study.values) for _,g in rhr.groupby('id')])
        ct=pd.read_csv(f'{data_dir}/computed_temperature.csv'); ct=ct[ct.type=='SKIN'].sort_values(['id','sleep_start_day_in_study']); report('nightly skin temperature (OBJECTIVE)',[regular_grid(g.nightly_temperature.values,g.sleep_start_day_in_study.values) for _,g in ct.groupby('id')])
    except FileNotFoundError: print('[objective wearable files missing; skipped]')
    out=pd.DataFrame(rows); sx=out.loc[out.measure.str.contains('PREDICTOR'),'spectral'].iloc[0]; sy=out.loc[out.measure.str.startswith(('fatigue','moodswing')),'spectral'].median(); proxy=np.sqrt(sx*sy)
    print(f'Identification note: sqrt(SSF_x*SSF_y)={proxy:.3f} is the f=0 sensitivity endpoint, NOT an estimate of reliability attenuation.')
    return out

def differential_prediction(df,rng):
    rows=[]
    for item,cls in CLASSIFICATION.items():
        per=paired_series(df,item)
        if len(per)<10: continue
        g=effect_g(pooled_within_corr(per)); pred='LEAK' if cls.startswith('DIRECTIONAL') else ('null' if cls.startswith('BALANCED') else '-'); obs='LEAKS' if g>=.10 else 'null'; rows.append(dict(item=item,cls=cls,group_g=g,predicted=pred,observed=obs))
    out=pd.DataFrame(rows); print(out.to_string(index=False)); return out

def attenuation_sensitivity(df,rng,n_boot=2000):
    sx=np.median([v for _,g in df.groupby('id') for v in [ssf_spectral(regular_grid(g.estrogen.values,g.day_in_study.values))] if np.isfinite(v)]); rows=[]
    for item,cls in CLASSIFICATION.items():
        if not (cls.startswith('BALANCED') or cls.startswith('DIRECTIONAL')): continue
        per=paired_series(df,item)
        if len(per)<10: continue
        sy=np.median([v for _,g in df.groupby('id') for v in [ssf_spectral(regular_grid(g[item].map(ORDINAL).values,g.day_in_study.values))] if np.isfinite(v)]); r_obs=pooled_within_corr(per)
        for f in (0.,.25,.50,.75):
            rx=sx+f*(1-sx); ry=sy+f*(1-sy); att=np.sqrt(rx*ry); boots=[]
            for _ in range(n_boot):
                idx=rng.integers(0,len(per),len(per)); pb=[per[i] for i in idx]; boots.append(effect_g(np.clip(pooled_within_corr(pb)/att,-.99,.99)))
            lo,hi=np.percentile(boots,[2.5,97.5]); rows.append(dict(item=item,cls=cls,f=f,ssf_predictor=sx,ssf_outcome=sy,implied_reliability_predictor=rx,implied_reliability_outcome=ry,attenuation=att,within_r_observed=r_obs,g_sensitivity=effect_g(np.clip(r_obs/att,-.99,.99)),ci95_low=lo,ci95_high=hi))
    out=pd.DataFrame(rows); print('No f is privileged by the data; attenuation is an identification sensitivity analysis.'); return out

def phase_locked(df):
    rows=[]
    for item in ('fatigue','moodswing','cramps','bloating','sorebreasts'):
        d=df.dropna(subset=['estrogen',item]).copy(); d['y']=d[item].map(ORDINAL); d=d.dropna(subset=['y']); d['yc']=d.y-d.groupby('id').y.transform('mean'); d['ec']=d.estrogen-d.groupby('id').estrogen.transform('mean'); r_e2=abs(corr(d.ec,d.yc)); ss_b=sum(len(g)*(g.mean()-d.yc.mean())**2 for _,g in d.groupby('phase').yc); ss_t=((d.yc-d.yc.mean())**2).sum(); eta2=ss_b/ss_t if ss_t>0 else np.nan; mens=d[d.phase.astype(str).str.contains('enstrual',case=False,na=False)].yc.mean(); rows.append(dict(item=item,r_e2_level=r_e2,phase_eta2=eta2,menstrual_mean=mens))
    out=pd.DataFrame(rows); print(out.to_string(index=False)); return out

def modelfree_bound(df,rng,item='fatigue',B=500):
    per=paired_series_with_days(df,item); r_obs=np.array([corr(x,y) for x,y,_ in per]); S=float(np.std(r_obs)); Sn=np.array([np.std([corr(phase_randomize_indexed(x,day,rng),y) for x,y,day in per]) for _ in range(B)]); p=(1+int(np.sum(Sn>=S)))/(B+1); print(f'{item}: SD(r_i)={S:.4f}, null median={np.median(Sn):.4f}, p={p:.3f}'); return dict(item=item,n=len(per),SD_obs=S,SD_null_median=float(np.median(Sn)),p=p)
def objective_vs_selfreport(df,data_dir,rng,B=300):
    rhr=pd.read_csv(f'{data_dir}/resting_heart_rate.csv')[['id','day_in_study','value']].rename(columns={'value':'rhr'}); ct=pd.read_csv(f'{data_dir}/computed_temperature.csv'); ct=ct[ct.type=='SKIN'][['id','sleep_start_day_in_study','nightly_temperature']].rename(columns={'sleep_start_day_in_study':'day_in_study','nightly_temperature':'temp'}); d=df.copy(); d['fatigue_n']=d.fatigue.map(ORDINAL); m=d[['id','day_in_study','estrogen','fatigue_n']].merge(rhr,on=['id','day_in_study'],how='left').merge(ct,on=['id','day_in_study'],how='left')
    complete={'fatigue_n':{pid:regular_grid(g.fatigue.map(ORDINAL).values,g.day_in_study.values) for pid,g in df.groupby('id')},'rhr':{pid:regular_grid(g.rhr.values,g.day_in_study.values) for pid,g in rhr.groupby('id')},'temp':{pid:regular_grid(g.temp.values,g.day_in_study.values) for pid,g in ct.groupby('id')}}; rows=[]
    for col,lab,kind in [('fatigue_n','fatigue','SELF-REPORT'),('rhr','resting heart rate','OBJECTIVE'),('temp','skin temperature','OBJECTIVE')]:
        per=[]; ss=[]
        for pid,g in m.groupby('id'):
            gg=g.dropna(subset=['estrogen',col]).sort_values('day_in_study')
            if len(gg)>=MIN_PAIRED:
                per.append((gg.estrogen.to_numpy(float),gg[col].to_numpy(float),gg.day_in_study.to_numpy(int))); y=complete[col].get(pid); v=ssf_spectral(y) if y is not None else np.nan
                if np.isfinite(v): ss.append(v)
        if len(per)<10: continue
        S=np.std([corr(x,y) for x,y,_ in per]); Sn=np.array([np.std([corr(phase_randomize_indexed(x,day,rng),y) for x,y,day in per]) for _ in range(B)]); p=(1+int(np.sum(Sn>=S)))/(B+1); rows.append(dict(outcome=lab,kind=kind,ssf=float(np.median(ss)) if ss else np.nan,SD_ri=S,null_SD=float(np.median(Sn)),p=p))
    out=pd.DataFrame(rows); print(out.to_string(index=False)); print('Objective comparisons are exploratory; multiplicity and progesterone confounding must travel with interpretation.'); return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--data-dir',required=True); ap.add_argument('--out-dir',default='./results'); ap.add_argument('--seed',type=int,default=11); ap.add_argument('--boot',type=int,default=2000); ap.add_argument('--surrogates',type=int,default=500); a=ap.parse_args(); os.makedirs(a.out_dir,exist_ok=True); rng=np.random.default_rng(a.seed); df=load(a.data_dir); describe(df)
    instrument_ssf(df,a.data_dir).to_csv(f'{a.out_dir}/table08_instrument_ssf.csv',index=False); differential_prediction(df,rng).to_csv(f'{a.out_dir}/table11_differential.csv',index=False); attenuation_sensitivity(df,rng,a.boot).to_csv(f'{a.out_dir}/table12_attenuation_sensitivity.csv',index=False); phase_locked(df).to_csv(f'{a.out_dir}/fig04_phase_locked.csv',index=False); pd.DataFrame([modelfree_bound(df,rng,'fatigue',a.surrogates),modelfree_bound(df,rng,'moodswing',a.surrogates)]).to_csv(f'{a.out_dir}/modelfree_bound.csv',index=False); objective_vs_selfreport(df,a.data_dir,rng,min(a.surrogates,500)).to_csv(f'{a.out_dir}/table14_objective.csv',index=False); print(f'[saved] {a.out_dir}/')
if __name__=='__main__': main()
