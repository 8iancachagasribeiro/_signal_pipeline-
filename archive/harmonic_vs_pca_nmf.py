#!/usr/bin/env python3
"""
harmonic_vs_pca_nmf.py — Does per-channel cyclic (harmonic) regression beat PCA and NMF
                         at recovering menstrual cycle phase?

THE PROPOSAL BEING TESTED
-------------------------
Fit each channel independently:      y = b0 + b1*cos(2*pi*t/T) + b2*sin(2*pi*t/T) + e
Per-channel phase:                   phi_j = atan2(b2, b1)
Per-channel amplitude:               A_j   = sqrt(b1^2 + b2^2)
Combine channels weighted by fit quality, then read the cycle phase directly.

This is harmonic regression; with irregular sampling it is the Lomb-Scargle
construction. It is established, not novel, which is an advantage: its properties
are known and it does not need to be validated from scratch.

THE CLAIM THAT MATTERS MOST
---------------------------
"Works with missing data."

That claim, if true, dissolves the blocker that killed the fusion analysis in
wearable_fusion.py. The FFT and PCA both assume complete, regularly spaced series.
mcPHASES wearable channels contain gaps, and requiring contiguous coverage of all six
channels leaves N = 1. Regression has no such requirement: it fits at whatever
timepoints exist.

WHAT THIS SCRIPT TESTS
----------------------
Four regimes crossed:
  shape : pure sinusoid  vs  realistic two-peak estradiol (ovulatory surge + luteal rise)
  gaps  : complete       vs  30% of days missing at random

against four methods: PCA, NMF (min-max scaled), harmonic regression with one
harmonic, and harmonic regression with two harmonics.

METRICS
-------
  |r| with the true latent   — does the recovered signal track the real cycle?
  phase error, in days       — the quantity the product actually needs

NOTE ON FAIRNESS
----------------
Under gaps, PCA and NMF are given complete-case rows (gap-free days only). They are
NOT given interpolated data, because interpolation inflates the smooth-signal fraction
by +0.07 to +0.10 at a 12% gap rate — see wearable_fusion.py, section 3. Dropping rows
is the honest comparison.

USAGE
-----
    python harmonic_vs_pca_nmf.py
"""
import warnings

import numpy as np
from sklearn.decomposition import NMF

warnings.filterwarnings("ignore")

CYCLE = 28.0
T_DAYS = 90
N_SUB = 60
SEED = 2026

BASE = np.array([62.0, 34.5, 33.8, 40.0, 800.0, 600.0])
AMPL = np.array([2.5, 0.35, 0.30, 4.0, 60.0, 45.0])
SIGNS = np.array([+1, +1, +1, -1, -1, +1])
SSF_TARGET = np.array([0.47, 0.34, 0.34, 0.40, 0.35, 0.35])


def latent_signal(t, shape):
    """The true underlying cycle."""
    if shape == "sinusoid":
        s = np.sin(2 * np.pi * t / CYCLE)
    else:
        # realistic estradiol: narrow ovulatory surge + broader luteal rise
        d = np.mod(t, CYCLE)
        s = (np.exp(-0.5 * ((d - 13) / 2.0) ** 2)
             + 0.55 * np.exp(-0.5 * ((d - 21) / 3.5) ** 2))
    return (s - s.mean()) / s.std()


def make_subject(rng, shape):
    t = np.arange(T_DAYS, dtype=float)
    lat = latent_signal(t, shape)
    X = np.empty((T_DAYS, 6))
    for j in range(6):
        sig = SIGNS[j] * AMPL[j] * lat
        s = SSF_TARGET[j]
        sd = np.sqrt(max(sig.var(), 1e-12) * (1 - s) / s)
        X[:, j] = BASE[j] + sig + rng.normal(0, sd, T_DAYS)
    return np.clip(X, 1e-9, None), lat, t


# ----------------------------------------------------------------- methods #
def pca_recover(X):
    sd = X.std(0)
    Z = (X - X.mean(0)) / np.where(sd < 1e-9, 1.0, sd)
    Zc = Z - Z.mean(0)
    U, S, Vt = np.linalg.svd(Zc, full_matrices=False)
    return Zc @ Vt.T[:, 0]


def nmf_recover(X, lat, k=2):
    lo, hi = X.min(0), X.max(0)
    A = np.clip((X - lo) / np.where(hi - lo < 1e-12, 1.0, hi - lo), 1e-9, None)
    W = NMF(n_components=k, init="nndsvda", max_iter=600,
            random_state=0, tol=1e-5).fit_transform(A)
    best_r, best = 0.0, W[:, 0]
    for j in range(k):
        if np.std(W[:, j]) < 1e-12:
            continue
        r = abs(np.corrcoef(W[:, j], lat)[0, 1])
        if r > best_r:
            best_r, best = r, W[:, j]
    return best


def harmonic_recover(X, t, n_harm=1):
    """Per-channel harmonic regression, combined by circular weighted mean of phase.

    Weights are R^2/(1-R^2) rather than amplitude*R^2: amplitude carries channel units
    (bpm vs degrees vs ms) and is not comparable across channels without standardising,
    whereas R^2/(1-R^2) is the signal-to-noise ratio of the fit and is unitless.
    """
    cols = [np.ones_like(t)]
    for h in range(1, n_harm + 1):
        cols += [np.cos(2 * np.pi * h * t / CYCLE), np.sin(2 * np.pi * h * t / CYCLE)]
    D = np.column_stack(cols)

    phases, weights = [], []
    for j in range(X.shape[1]):
        y = X[:, j]
        beta, *_ = np.linalg.lstsq(D, y, rcond=None)
        resid = y - D @ beta
        ss_tot = ((y - y.mean()) ** 2).sum()
        r2 = 1.0 - (resid ** 2).sum() / ss_tot if ss_tot > 0 else 0.0
        r2 = float(np.clip(r2, 0.0, 0.999))
        b1, b2 = beta[1], beta[2]          # first harmonic
        phases.append(np.arctan2(b2, b1))
        weights.append(r2 / (1.0 - r2))    # SNR of the fit

    phases, weights = np.array(phases), np.array(weights)
    if weights.sum() <= 0:
        return np.zeros_like(t)
    # circular weighted mean, respecting that a sign flip is a phase shift of pi
    ref = phases[np.argmax(weights)]
    aligned = phases + np.pi * (np.abs(np.angle(np.exp(1j * (phases - ref)))) > np.pi / 2)
    C = np.sum(weights * np.cos(aligned))
    S = np.sum(weights * np.sin(aligned))
    phi = np.arctan2(S, C)
    return np.cos(2 * np.pi * t / CYCLE - phi)


# ----------------------------------------------------------------- scoring #
def phase_error_days(rec, lat, t):
    """Circular distance between recovered and true phase, expressed in days.

    CORRECTION (v2). An earlier version of this function estimated phase with an FFT.
    That was wrong under the gapped regimes: the FFT assumes regular spacing, and
    applying it to complete-case rows treats non-consecutive days as consecutive. It is
    the same error this project's companion manuscript documents in section 7.8, and it
    invalidated every phase figure in the rows with missing data.

    The estimator below fits a single harmonic at the cycle frequency BY LEAST SQUARES
    ON THE ACTUAL TIMEPOINTS. This is the Lomb-Scargle construction, and it is valid for
    arbitrary sampling, regular or not. It is applied identically to every method, so the
    comparison stays fair.
    """
    def phase_of(y, t):
        D = np.column_stack([np.ones_like(t),
                             np.cos(2 * np.pi * t / CYCLE),
                             np.sin(2 * np.pi * t / CYCLE)])
        beta, *_ = np.linalg.lstsq(D, y - np.mean(y), rcond=None)
        return np.arctan2(beta[2], beta[1])

    d = np.angle(np.exp(1j * (phase_of(rec, t) - phase_of(lat, t))))
    return abs(d) / (2 * np.pi) * CYCLE


def run(shape, gap_frac):
    rng = np.random.default_rng(SEED)
    out = {m: {"r": [], "ph": []} for m in
           ("PCA", "NMF", "harmonic k=1", "harmonic k=2")}

    for _ in range(N_SUB):
        X, lat, t = make_subject(rng, shape)

        if gap_frac > 0:
            keep = rng.random(T_DAYS) > gap_frac
            if keep.sum() < 30:
                continue
            Xg, tg, latg = X[keep], t[keep], lat[keep]
        else:
            Xg, tg, latg = X, t, lat

        # PCA and NMF get complete-case rows only — no interpolation
        for name, rec in (("PCA", pca_recover(Xg)),
                          ("NMF", nmf_recover(Xg, latg)),
                          ("harmonic k=1", harmonic_recover(Xg, tg, 1)),
                          ("harmonic k=2", harmonic_recover(Xg, tg, 2))):
            if np.std(rec) < 1e-12:
                continue
            out[name]["r"].append(abs(np.corrcoef(rec, latg)[0, 1]))
            out[name]["ph"].append(phase_error_days(rec, latg, tg))

    tag = f"{shape.upper()} shape | {int(gap_frac*100)}% days missing"
    print(f"\n{tag}")
    print("-" * 62)
    print(f"{'method':>14} {'|r| with truth':>16} {'phase error (days)':>22}")
    for m, v in out.items():
        if not v["r"]:
            print(f"{m:>14} {'FAILED':>16}")
            continue
        print(f"{m:>14} {np.nanmedian(v['r']):>16.3f} {np.nanmedian(v['ph']):>22.2f}")


def main():
    print("=" * 62)
    print("HARMONIC REGRESSION vs PCA vs NMF — cycle phase recovery")
    print("=" * 62)
    print(f"  {N_SUB} subjects, {T_DAYS} days, 6 mixed-sign non-negative channels")
    print("  noise calibrated to the smooth-signal fractions measured in mcPHASES")

    for shape in ("sinusoid", "two-peak"):
        for gap in (0.0, 0.30):
            run(shape, gap)

    print("\n" + "=" * 62)
    print("  Phase error is the quantity the product needs. Below ~2 days is")
    print("  usable for phase-conditioned prediction; above ~4 days is not.")


if __name__ == "__main__":
    main()
