#!/usr/bin/env python3
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path


def run(cmd):
    print('[run]', ' '.join(map(str,cmd)), flush=True)
    subprocess.run(cmd, check=True)


def main():
    ap=argparse.ArgumentParser(description='Reproduce the BRM methodological expansion.')
    ap.add_argument('--root',default='.')
    ap.add_argument('--reps',type=int,default=1000)
    ap.add_argument('--ssf-workers',type=int,default=4)
    ap.add_argument('--surrogate-workers',type=int,default=3)
    a=ap.parse_args()
    root=Path(a.root).resolve(); results=root/'results'; results.mkdir(parents=True,exist_ok=True)
    py=sys.executable
    run([py,str(root/'reference_model.py')])
    run([py,str(root/'brm_robustness.py'),'--reps',str(a.reps),'--out-dir',str(results)])
    run([py,str(root/'brm_sampling_design.py'),'--reps',str(a.reps),'--out-dir',str(results)])
    run([py,str(root/'brm_ssf_benchmark.py'),'--reps',str(a.reps),'--workers',str(a.ssf_workers),'--out-dir',str(results)])
    run([py,str(root/'brm_surrogate_power.py'),'--reps',str(a.reps),'--surrogates','199','--workers',str(a.surrogate_workers),'--out-dir',str(results)])
    run([py,str(root/'make_brm_outputs.py'),'--root',str(root),'--results-dir','results','--figures-dir','figures','--tables-dir','tables'])
    if a.reps == 1000:
        run([py,str(root/'validate_brm_outputs.py'),'--root',str(root)])
    else:
        print('[skip] canonical validator requires --reps 1000', flush=True)

if __name__=='__main__': main()
