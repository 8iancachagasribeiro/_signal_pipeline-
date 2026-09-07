#!/usr/bin/env python3
"""Independent clinical-actigraphy SSF replication.

HYPERAKTIV labels are read from patient_info.csv because activity_data mixes ADHD
participants and controls. Unresolved IDs are excluded, never guessed. Gap-free hourly
segments are used for FFT-based SSF, and outputs are saved so the final N is derived
from the current run rather than copied from an older manuscript.
"""
import argparse,glob,os,re,warnings
import numpy as np,pandas as pd
from ssf_estimators import ssf_spectral
warnings.filterwarnings('ignore')
SOURCES=[('depresjon/data/condition/*.csv','major depression'),('depresjon/data/control/*.csv','controls'),('psykose/patient/*.csv','schizophrenia'),('psykose/control/*.csv','controls')]

def load_series(path):
    d=None
    for sep in (',',';'):
        try:
            q=pd.read_csv(path,sep=sep)
            if q.shape[1]>=2: d=q; break
        except Exception: pass
    if d is None: return None
    acts=[c for c in d.columns if 'activ' in c.lower()]; times=[c for c in d.columns if c.lower() in ('timestamp','date','time')]
    if not acts or not times: return None
    t=pd.to_datetime(d[times[0]],errors='coerce',format='mixed'); y=pd.to_numeric(d[acts[0]],errors='coerce')
    return pd.Series(y.values,index=t).dropna()

def longest_contiguous(s,step_hours=1.):
    if len(s)<2: return s
    gaps=pd.Series(s.index).diff().dt.total_seconds().fillna(step_hours*3600).values/3600.; best_len=best_start=start=0
    for k in range(1,len(gaps)+1):
        if k==len(gaps) or abs(gaps[k]-step_hours)>1e-6:
            if k-start>best_len: best_len,best_start=k-start,start
            start=k
    return s.iloc[best_start:best_start+best_len]

def _process(path,label,seen,rows,min_hours):
    s=load_series(path)
    if s is None or len(s)<min_hours: return False
    s=s.resample('1h').mean().dropna(); seg=longest_contiguous(s)
    if len(seg)<min_hours: return False
    key=(round(float(seg.mean()),4),len(seg))
    if key in seen: return False
    seen.add(key); v=ssf_spectral(seg.values)
    if not np.isfinite(v): return False
    rows.append(dict(group=label,hours=len(seg),ssf=v,series=seg.values)); return True

def _find_metadata(root,explicit=None):
    cand=[]
    if explicit: cand.append(explicit if os.path.isabs(explicit) else os.path.join(root,explicit))
    cand += [os.path.join(root,'patient_info.csv'),os.path.join(root,'hyperaktiv','patient_info.csv')]
    for p in cand:
        if os.path.isfile(p): return p
    hits=glob.glob(os.path.join(root,'**','patient_info.csv'),recursive=True); return hits[0] if hits else None

def _id(path):
    m=re.search(r'(\d+)',os.path.basename(path)); return int(m.group(1)) if m else None

def _adhd_ids(root,metadata=None,adhd_column=None,adhd_positive=None):
    p=_find_metadata(root,metadata)
    if p is None: raise FileNotFoundError('HYPERAKTIV patient_info.csv not found; refusing to infer diagnosis from folder membership')
    m=pd.read_csv(p,sep=None,engine='python'); m.columns=[str(c).strip() for c in m.columns]
    idcol=next((c for c in ('ID','id','Id','subject','SUBJECT','SubjectID') if c in m.columns),None)
    if idcol is None: raise KeyError(f'No ID column in {p}')
    if adhd_column is None:
        if 'ADHD' in m.columns: adhd_column='ADHD'; kind='num'
        else:
            adhd_column=next((c for c in ('DIAGNOSIS','diagnosis','GROUP','group','label','LABEL') if c in m.columns),None); kind='str'
            if adhd_column is None: raise KeyError(f'No ADHD label column in {p}')
    else: kind='str'
    target=str(adhd_positive).strip().upper() if adhd_positive is not None else None
    ids=set()
    for _,r in m.iterrows():
        try:
            positive=(str(r[adhd_column]).strip().upper()==target) if target is not None else ((float(r[adhd_column])==1.) if kind=='num' else str(r[adhd_column]).strip().upper()=='ADHD')
            if positive: ids.add(int(float(r[idcol])))
        except Exception: pass
    if not ids: raise ValueError('Parsed 0 ADHD subjects; refusing to proceed')
    return ids,os.path.basename(p)

def collect(root,min_hours=48,metadata=None,adhd_column=None,adhd_positive=None):
    rows=[]; seen=set()
    for pattern,label in SOURCES:
        for p in sorted(glob.glob(os.path.join(root,pattern))): _process(p,label,seen,rows,min_hours)
    files=sorted(glob.glob(os.path.join(root,'activity_data/*.csv')))
    if files:
        ids,source=_adhd_ids(root,metadata,adhd_column,adhd_positive); kept=excluded=unresolved=0
        for p in files:
            sid=_id(p)
            if sid is None: unresolved+=1; continue
            if sid not in ids: excluded+=1; continue
            kept += int(_process(p,'ADHD',seen,rows,min_hours))
        print(f'[HYPERAKTIV] labels={source}; files={len(files)} ADHD kept={kept} non-ADHD excluded={excluded} unresolved={unresolved}')
    return rows

def table15(rows):
    df=pd.DataFrame([{k:r[k] for k in ('group','hours','ssf')} for r in rows])
    print(df.groupby('group').agg(n=('ssf','size'),hours_median=('hours','median'),ssf_median=('ssf','median')).to_string())
    print(f'\nTOTAL n={len(df)}; median SSF={df.ssf.median():.3f}')
    print('Cross-domain contrast is descriptive: construct, device, aggregation, missingness and cadence differ; higher circadian SSF does not prove equal device precision across domains.')
    return df

def table16(rows,min_cycles=16):
    long=[r['series'] for r in rows if len(r['series'])>=min_cycles*24]; out=[]
    for nc in (3,5,8,12,16):
        vals=np.asarray([ssf_spectral(s[:nc*24]) for s in long],float); vals=vals[np.isfinite(vals)]
        out.append(dict(cycles=nc,samples=nc*24,ssf_median=float(np.median(vals)) if len(vals) else np.nan,q25=float(np.percentile(vals,25)) if len(vals) else np.nan,q75=float(np.percentile(vals,75)) if len(vals) else np.nan,n=len(vals)))
    df=pd.DataFrame(out); print(df.to_string(index=False)); return df

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--data-dir',required=True); ap.add_argument('--min-hours',type=int,default=48); ap.add_argument('--hyperaktiv-metadata',default=None); ap.add_argument('--adhd-column',default=None); ap.add_argument('--adhd-positive',default=None); ap.add_argument('--out-dir',default='./results'); a=ap.parse_args()
    rows=collect(a.data_dir,a.min_hours,a.hyperaktiv_metadata,a.adhd_column,a.adhd_positive)
    if not rows: raise SystemExit('no usable series found')
    os.makedirs(a.out_dir,exist_ok=True); t15=table15(rows); t16=table16(rows); t15.to_csv(f'{a.out_dir}/table15_actigraphy_ssf.csv',index=False); t16.to_csv(f'{a.out_dir}/table16_cycle_count_robustness.csv',index=False); print(f'[saved] {a.out_dir}/')
if __name__=='__main__': main()
