#!/usr/bin/env python3
"""Utilities for SSF-calibrated simulation and power sensitivity analyses.

Smooth-signal fraction (SSF) is a variance fraction, not psychometric reliability.
Power based only on SSF is not point-identifiable because SSF does not reveal what
share of smooth outcome variance is actually driven by the predictor of interest.
"""
from __future__ import annotations
import numpy as np


def _standardize(x):
    x=np.asarray(x,float); x=x-np.nanmean(x); sd=np.nanstd(x)
    return np.zeros_like(x) if (not np.isfinite(sd) or sd<1e-12) else x/sd


def add_white_noise_for_ssf(smooth,target_ssf,rng):
    if not (0<target_ssf<=1): raise ValueError('target_ssf must lie in (0, 1]')
    s=_standardize(smooth)
    if target_ssf==1: return s
    return s+rng.normal(0,np.sqrt((1-target_ssf)/target_ssf),len(s))


def smooth_nuisance(n,rng,max_cycles_per_sample=0.12):
    n=int(n)
    if n<4: return _standardize(rng.normal(size=n))
    f=np.fft.rfftfreq(n,d=1.0); spec=np.zeros(len(f),dtype=complex)
    allowed=np.where((f>0)&(f<=max_cycles_per_sample))[0]
    if len(allowed)==0: return _standardize(np.linspace(-1,1,n))
    amp=1.0/np.sqrt(np.maximum(f[allowed],1.0/n)); phase=rng.uniform(0,2*np.pi,len(allowed))
    spec[allowed]=amp*np.exp(1j*phase); z=np.fft.irfft(spec,n=n)
    z=z+rng.normal(0,.15)*np.linspace(-1,1,n)
    return _standardize(z)


def outcome_from_mechanism(mechanistic_smooth,target_ssf,coupled_fraction,rng):
    q=float(coupled_fraction)
    if not (0<=q<=1): raise ValueError('coupled_fraction must lie in [0, 1]')
    m=_standardize(mechanistic_smooth); z=smooth_nuisance(len(m),rng)
    smooth=_standardize(np.sqrt(q)*m+np.sqrt(1-q)*z)
    return add_white_noise_for_ssf(smooth,target_ssf,rng)


def add_white_noise_for_estimated_ssf(smooth,target_ssf,rng,estimator,tol=0.003,max_iter=32):
    """Scale one fixed white-noise draw until the adopted estimator hits target SSF."""
    if not (0<target_ssf<=1): raise ValueError('target_ssf must lie in (0, 1]')
    s=_standardize(smooth)
    if target_ssf==1: return s
    z=_standardize(rng.normal(size=len(s))); lo,hi=0.0,10.0
    best=add_white_noise_for_ssf(s,target_ssf,rng); best_err=np.inf
    for _ in range(max_iter):
        mid=(lo+hi)/2; y=s+mid*z; est=estimator(y)
        if not np.isfinite(est): break
        err=abs(est-target_ssf)
        if err<best_err: best,best_err=y.copy(),err
        if err<=tol: break
        if est>target_ssf: lo=mid
        else: hi=mid
    return best


def summarize_actual_ssf(series_list,estimator):
    vals=np.asarray([estimator(np.asarray(x,float)) for x in series_list],float)
    vals=vals[np.isfinite(vals)]
    return float(np.median(vals)) if len(vals) else np.nan
