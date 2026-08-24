#!/usr/bin/env python3
"""Static + smoke-test audit for the Article 1 reproducibility package."""
from __future__ import annotations
import ast, glob, os
import numpy as np

ROOT=os.path.dirname(os.path.abspath(__file__))
FAIL=[]


def check(cond,msg):
    (print('[PASS]',msg) if cond else (FAIL.append(msg),print('[FAIL]',msg)))


def static_checks():
    py=glob.glob(os.path.join(ROOT,'*.py'))
    bad=[]
    for f in py:
        txt=open(f,encoding='utf-8').read()
        try: ast.parse(txt)
        except SyntaxError as e: bad.append(f'{os.path.basename(f)}:{e.lineno}')
    check(not bad,'all Python files parse')
    forbidden=[]
    for f in [x for x in py+glob.glob(os.path.join(ROOT,'*.sh')) if os.path.basename(x) != 'audit_consistency.py']:
        txt=open(f,encoding='utf-8',errors='ignore').read()
        for token in ('/home/claude','/mnt/user-data/outputs'):
            if token in txt: forbidden.append((os.path.basename(f),token))
    check(not forbidden,'no environment-specific output paths remain')
    mf=open(os.path.join(ROOT,'make_figures.py'),encoding='utf-8').read()
    check(all(f'def fig{i}(' in mf for i in range(1,6)),'figure script implements all five current figures')


def smoke_ssf():
    from registered_test_power import ssf_diagnostic
    rng=np.random.default_rng(123)
    sx,sy=ssf_diagnostic(rng,q=.50)
    check(abs(sx-.469)<.015 and abs(sy-.323)<.015,
          f'SSF-calibrated simulator hits targets (got {sx:.3f}/{sy:.3f})')


def smoke_mechanism():
    import h4_frontier as H
    rng=np.random.default_rng(0); gr=[]; indiv=[]; pos=[]
    for _ in range(200):
        df=H.simulate_study(rng,H.N_SUBJ,28,2,.73,.05)
        r=H.within_person_r(df,'y_true'); rr=H.group_effect(df,'y_true')
        gr.append(rr); indiv.append(np.median(np.abs(r))); pos.append(np.mean(r>0))
    rg=float(np.mean(gr)); g=abs(2*rg/np.sqrt(max(1-rg*rg,1e-12)))
    ri=float(np.mean(indiv)); gi=2*ri/np.sqrt(max(1-ri*ri,1e-12))
    check(.04 <= g <= .12 and .25 <= gi <= .40 and .35 <= np.mean(pos) <= .60,
          f'core masking calibration remains in preregistered neighborhood (group |g|={g:.3f}, individual |g|={gi:.3f})')


def smoke_varying_lrt():
    from fastlrt import lrt_random_slope_fast,lrt_random_slope_varying_x
    rng=np.random.default_rng(1); N=10;n=18;x=np.linspace(-1,1,n); sl=rng.normal(0,.5,N)
    Y=np.array([sl[i]*x+rng.normal(0,.5,n) for i in range(N)]); X=np.tile(x,(N,1))
    p1=lrt_random_slope_fast(Y,x); p2=lrt_random_slope_varying_x(Y,X)
    check(np.isfinite(p1) and np.isfinite(p2) and abs(p1-p2)<1e-5,
          'varying-predictor LRT reproduces shared-predictor LRT when X is identical')


def static_semantics():
    txt=open(os.path.join(ROOT,'registered_test_power.py'),encoding='utf-8').read()
    check('power is a sensitivity surface' in txt and 'coupled_fraction' in txt,
          'power script labels SSF analysis as sensitivity, not empirical point power')
    txt2=open(os.path.join(ROOT,'budget_allocation.py'),encoding='utf-8').read()
    check('No single row is' in txt2,'budget script rejects a single empirical-power interpretation')


def main():
    static_checks(); smoke_ssf(); smoke_mechanism(); smoke_varying_lrt(); static_semantics()
    print('\nOVERALL:', 'PASS' if not FAIL else 'FAIL')
    if FAIL:
        for x in FAIL: print(' -',x)
        raise SystemExit(1)

if __name__=='__main__': main()
