#!/usr/bin/env python3
"""
nmf_vs_pca.py — Does non-negative matrix factorisation recover a cyclic latent
                better than PCA on non-negative physiological channels?

WHY THIS TEST EXISTS
--------------------
 observed that wearable signals are non-negative and suggested
NMF may therefore suit the latent-space extraction better than PCA.

The observation is correct about the RAW signals. But two things need testing:

  (a) The PCA in wearable_fusion.py consumes Z-SCORED data, which is NOT non-negative.
      So the premise does not describe what was actually run.

  (b) More importantly: the target is a cyclic OSCILLATION about a baseline. NMF
      factorises X ~ WH with W,H >= 0, i.e. as additive combinations of non-negative
      parts. Representing a downswing requires a negative contribution, which NMF
      forbids. The prediction is that NMF's leading component is dominated by the
      BASELINE LEVEL, not the oscillation.

  (c) And decisively: real channels couple to the cycle with DIFFERENT SIGNS
      (skin temperature rises in the luteal phase; other measures fall). A single
      non-negative component cannot represent a sign flip across channels.

This script tests (b) and (c) against ground truth, which mcPHASES cannot provide:
requiring contiguous gap-free coverage of all six channels there leaves N = 1.

DESIGN
------
Six synthetic channels with realistic baselines and amplitudes, all non-negative,
each carrying the SAME latent cycle with a channel-specific amplitude and sign.
Noise is calibrated to the smooth-signal fractions actually measured in mcPHASES
(0.32 to 0.47).

Two regimes:
  ALIGNED  — every channel couples to the cycle with the same sign (NMF's best case)
  MIXED    — signs differ across channels (the physiologically realistic case)

METRICS
-------
  |r| with the true latent  : does the component track the real cycle?
  SSF of the component      : is the component itself cycle-recoverable?

USAGE
-----
    python nmf_vs_pca.py
"""
import warnings

import numpy as np
from sklearn.decomposition import NMF

warnings.filterwarnings("ignore")

LN2 = np.log(2.0)
CYCLE = 28.0
T = 90                      # days, matching the mcPHASES median
N_SUB = 60
SEED = 2026

# realistic baselines / amplitudes for the six mcPHASES wearable channels
BASE = np.array([62.0, 34.5, 33.8, 40.0, 800.0, 600.0])   # RHR, skin T, wrist T, RMSSD, LF, HF
AMPL = np.array([2.5, 0.35, 0.30, 4.0, 60.0, 45.0])       # cycle-linked excursion
SIGN_MIXED = np.array([+1, +1, +1, -1, -1, +1])
SIGN_ALIGNED = np.ones(6)
SSF_TARGET = np.array([0.47, 0.34, 0.34, 0.40, 0.35, 0.35])


def ssf_spectral(y, f_cut=0.25):
    """Adopted estimator: spectral separation with exponential-median bias correction."""
    y = np.asarray(y, float)
    y = y[np.isfinite(y)]
    n = len(y)
    if n < 25 or np.std(y) < 1e-12:
        return np.nan
    y = y - y.mean()
    P = (np.abs(np.fft.rfft(y)) ** 2) / n
    f = np.fft.rfftfreq(n, d=1.0)
    P = P[1:]
    f = f[1:]
    hi = f > f_cut
    if hi.sum() < 4:
        return np.nan
    return float(np.clip(1.0 - (np.median(P[hi]) / LN2) * len(P) / P.sum(), 0.0, 1.0))


def make_subject(rng, signs):
    """Six non-negative channels sharing one latent cycle."""
    t = np.arange(T, dtype=float)
    latent = np.sin(2 * np.pi * t / CYCLE)          # zero-mean, unit-ish amplitude
    X = np.empty((T, 6))
    for j in range(6):
        sig = signs[j] * AMPL[j] * latent
        # noise sd such that the smooth fraction of the channel matches SSF_TARGET
        s = SSF_TARGET[j]
        sd = np.sqrt(max(sig.var(), 1e-12) * (1 - s) / s)
        X[:, j] = BASE[j] + sig + rng.normal(0, sd, T)
    return np.clip(X, 1e-9, None), latent      # clipped: channels stay non-negative


def pca_component(X):
    """z-score, then leading principal component (what wearable_fusion.py actually does)."""
    sd = X.std(0)
    Z = (X - X.mean(0)) / np.where(sd < 1e-9, 1.0, sd)
    Zc = Z - Z.mean(0)
    U, S, Vt = np.linalg.svd(Zc, full_matrices=False)
    return Zc @ Vt.T[:, 0]


def nmf_components(X, k=2, scale="raw"):
    """NMF on non-negative input. Returns all k component time courses."""
    if scale == "minmax":
        lo, hi = X.min(0), X.max(0)
        A = (X - lo) / np.where(hi - lo < 1e-12, 1.0, hi - lo)
        A = np.clip(A, 1e-9, None)
    else:
        A = X
    m = NMF(n_components=k, init="nndsvda", max_iter=600, random_state=0, tol=1e-5)
    W = m.fit_transform(A)
    return W


def best_by_truth(comps, latent):
    """Pick the component that best tracks the true latent (gives NMF its best shot)."""
    best_r, best_c = 0.0, comps[:, 0]
    for j in range(comps.shape[1]):
        c = comps[:, j]
        if np.std(c) < 1e-12:
            continue
        r = abs(np.corrcoef(c, latent)[0, 1])
        if r > best_r:
            best_r, best_c = r, c
    return best_r, best_c


def run(signs, label):
    rng = np.random.default_rng(SEED)
    rows = {k: {"r": [], "ssf": []} for k in ("PCA", "NMF raw", "NMF minmax")}
    for _ in range(N_SUB):
        X, latent = make_subject(rng, signs)

        p = pca_component(X)
        rows["PCA"]["r"].append(abs(np.corrcoef(p, latent)[0, 1]))
        rows["PCA"]["ssf"].append(ssf_spectral(p))

        for tag, sc in (("NMF raw", "raw"), ("NMF minmax", "minmax")):
            W = nmf_components(X, k=2, scale=sc)
            r, c = best_by_truth(W, latent)
            rows[tag]["r"].append(r)
            rows[tag]["ssf"].append(ssf_spectral(c))

    print(f"\n{label}")
    print("-" * 58)
    print(f"{'method':>12} {'|r| with truth':>16} {'SSF of component':>19}")
    for k, v in rows.items():
        r = np.nanmedian(v["r"])
        s = np.nanmedian(v["ssf"])
        print(f"{k:>12} {r:>16.3f} {s:>19.3f}")
    return rows


def main():
    print("=" * 58)
    print("CAN NMF RECOVER A CYCLIC LATENT BETTER THAN PCA?")
    print("=" * 58)
    print(f"  {N_SUB} subjects, {T} days, 28-day cycle, 6 non-negative channels")
    print(f"  channel SSF calibrated to mcPHASES: {list(SSF_TARGET)}")

    run(SIGN_ALIGNED, "REGIME 1 — ALIGNED signs (NMF's best case)")
    run(SIGN_MIXED, "REGIME 2 — MIXED signs (physiologically realistic)")

    print("\n" + "=" * 58)
    print("WHAT TO READ")
    print("=" * 58)
    print("  |r| with truth near 1.0  => the component IS the cycle")
    print("  |r| near 0.0             => the component is something else (usually the")
    print("                              baseline level, which carries no cycle)")
    print("\n  If NMF loses in the MIXED regime, the reason is structural and not")
    print("  fixable by tuning: a single non-negative component cannot represent")
    print("  channels that move in opposite directions.")


if __name__ == "__main__":
    main()
