#!/usr/bin/env python3
"""
ssf_estimators.py — Smooth Signal Fraction (SSF) estimators, with validation.

WHAT THIS IS
------------
The SSF is the fraction of a series' observed variance carried by a temporally
STRUCTURED (smooth) component. For detecting a CYCLE-LOCKED coupling this is the
quantity that matters: no white component -- of any origin (measurement error,
genuine white state fluctuation, or ordinal quantisation) -- can carry a cycle-locked
signal. It is NOT classical psychometric reliability and must not be reported as such.

THREE ESTIMATORS
----------------
1. AR(1) closed form:      R = rho(1)^2 / rho(2)
   Assumes the true signal is AR(1). It is NOT: a cycle-locked symptom (e.g. menstrual
   cramps) is a periodic, sharply-peaked process. This estimator is BIASED and can
   return values > 1 (impossible for a variance fraction). It produced 1.063 for cramps
   in the v2 manuscript and was replaced.

2. ACF linear extrapolation to lag 0:   R = 2*rho(1) - rho(2)
   Better, but still biased: the autocorrelation of ANY differentiable process has zero
   derivative at lag 0 (it is even, with a maximum there), so rho(k) ~ 1 - a*k^2 --
   QUADRATIC, not linear.

3. SPECTRAL (adopted):
   White noise has a FLAT power spectrum; a cycle-locked signal (period ~28 samples) has
   no power at high frequency. The high-frequency plateau of the periodogram IS the noise
   floor. Makes NO assumption about the signal's shape.
   BIAS CORRECTION: periodogram ordinates of white noise are EXPONENTIALLY distributed,
   and the median of an exponential is ln(2) ~ 0.693 times its mean. Using the raw median
   underestimates the noise floor by ~30%. We divide by ln(2).

GAP HANDLING (added after audit; the reason the two domains are now symmetric)
------------------------------------------------------------------------------
The periodogram assumes REGULAR spacing. A gapped series must never be silently
compacted (treating non-consecutive samples as consecutive corrupts the spectrum and
biases the SSF). Two kinds of gap exist, and they need different handling:

  * NaN-marked gaps: a missing sample at a KNOWN grid position (value is NaN). The
    estimators see these directly and apply the `on_gap` policy below.
  * MISSING-ROW gaps: an absent row (e.g. a day with no record at all). A bare array
    CANNOT reveal these -- there is no NaN and no time index. The CALLER must expose
    them first with regular_grid(values, day_index), which reindexes onto the contiguous
    integer grid and inserts NaN at the holes. Only then is the `on_gap` policy able to
    act. If you pass an already-compacted array (missing rows simply absent), the
    estimator has no way to know, and the old bias returns. USE regular_grid().

`on_gap` policy (all three estimators):
    "longest" (default) -> compute on the longest contiguous gap-free run. Correct and
                           non-destructive; matches the actigraphy pipeline's
                           longest_contiguous() choice, so both domains use one rule.
    "error"             -> raise on any gap. Use in strict/audit runs to surface silent
                           gaps rather than absorb them.
    "compact"           -> LEGACY, UNSAFE: drop NaN and treat survivors as consecutive.
                           Reproduces the pre-audit behaviour; kept ONLY for explicit
                           before/after comparison (see validate_gaps()).

SCALE INVARIANCE
----------------
f_cut is in CYCLES PER SAMPLE, not cycles/day: the FFT uses rfftfreq(n, d=1.0). f_cut =
0.25 means "periods shorter than 4 SAMPLES are the noise band", whatever the sampling
interval is. On daily sampling that is < 4 days (the 28-day cycle, f = 0.036, is far
below); on hourly sampling it is < 4 hours (the 24-hour rhythm, f = 0.042, equally
below). This scale invariance is exactly what licenses the day<->hour transfer in the
cross-domain replication.

VALIDATION RESULT (validate(); three signal shapes x four noise levels, no gaps):
    estimator          mean |bias|   max |bias|
    AR(1)                0.077         0.157
    ACF-linear           0.036         0.082
    SPECTRAL (adopted)   0.028         0.082
validate_gaps() additionally exercises the gapped path the table above never touched.
"""
import warnings

import numpy as np

LN2 = np.log(2.0)


# ------------------------------- gap handling ------------------------------- #
def longest_finite_run(y):
    """Longest contiguous run of finite values in a REGULAR-grid array.

    NaN marks a missing sample at a known grid position. We do NOT compact across gaps
    (that would treat non-consecutive samples as consecutive and corrupt the
    periodogram); we return the longest gap-free stretch instead. No interpolation:
    filling gaps injects artificial smoothness and biases the SSF upward.
    """
    y = np.asarray(y, float)
    finite = np.isfinite(y)
    n = len(finite)
    best_len, best_start = 0, 0
    start = 0
    while start < n:
        if not finite[start]:
            start += 1
            continue
        end = start
        while end < n and finite[end]:
            end += 1
        if end - start > best_len:
            best_len, best_start = end - start, start
        start = end
    return y[best_start:best_start + best_len]


def regular_grid(values, index):
    """Place `values` on the contiguous integer grid implied by `index` (day_in_study,
    hour, ...), inserting NaN at missing positions. This is how a caller exposes
    MISSING-ROW gaps (absent rows, not NaN) to the estimator, which otherwise cannot
    see them from a bare array. Duplicate indices collapse to the first occurrence.
    Pass the result to any estimator below.
    """
    idx = np.asarray(index)
    val = np.asarray(values, float)
    if len(idx) != len(val):
        raise ValueError("values and index must have equal length")
    if len(idx) == 0:
        return val
    order = np.argsort(idx, kind="mergesort")
    idx, val = idx[order], val[order]
    keep = np.concatenate(([True], np.diff(idx) != 0))     # drop duplicate positions
    idx, val = idx[keep], val[keep]
    lo, hi = int(idx.min()), int(idx.max())
    grid = np.full(hi - lo + 1, np.nan)
    grid[(idx - lo).astype(int)] = val
    return grid


def _prep(y, on_gap):
    """Shared gap policy for all three estimators. See module docstring."""
    y = np.asarray(y, float)
    if np.isfinite(y).all():
        return y                                            # regular, gap-free: fast path
    if on_gap == "longest":
        return longest_finite_run(y)
    if on_gap == "error":
        raise ValueError(
            "series contains gaps (non-finite values). Pass a regular gap-free series, "
            "or regular_grid()+on_gap='longest'. Refusing to silently compact."
        )
    if on_gap == "compact":
        return y[np.isfinite(y)]                            # LEGACY, UNSAFE
    raise ValueError(f"unknown on_gap policy: {on_gap!r}")


# ------------------------------- estimators ------------------------------- #
def ssf_ar1(y, on_gap="longest"):
    """AR(1) closed form. BIASED for cyclic processes; can exceed 1. Reported for
    comparison only -- do not use for inference."""
    y = _prep(y, on_gap)
    if len(y) < 25 or np.std(y) < 1e-12:
        return np.nan
    y = y - y.mean(); n = len(y); v = np.dot(y, y) / n
    r1 = np.dot(y[:-1], y[1:]) / n / v
    r2 = np.dot(y[:-2], y[2:]) / n / v
    return r1 ** 2 / r2 if r2 > 0.02 else np.nan


def ssf_acf_linear(y, on_gap="longest"):
    """Linear extrapolation of the ACF to lag 0. Biased (true ACF is quadratic near 0)."""
    y = _prep(y, on_gap)
    if len(y) < 25 or np.std(y) < 1e-12:
        return np.nan
    y = y - y.mean(); n = len(y); v = np.dot(y, y) / n
    r1 = np.dot(y[:-1], y[1:]) / n / v
    r2 = np.dot(y[:-2], y[2:]) / n / v
    return 2 * r1 - r2


def ssf_spectral(y, f_cut=0.25, on_gap="longest"):
    """ADOPTED ESTIMATOR. Spectral separation with exponential-median bias correction.

    f_cut is in CYCLES PER SAMPLE (rfftfreq with d=1.0). 0.25 => periods shorter than
    4 SAMPLES are the noise band. Scale-invariant across daily/hourly sampling.
    `on_gap` governs gap handling; see module docstring. Callers with missing-ROW gaps
    must pass regular_grid(values, index) first, or the gaps stay invisible here.
    """
    y = _prep(y, on_gap)
    n = len(y)
    if n < 25 or np.std(y) < 1e-12:
        return np.nan
    y = y - y.mean()
    P = (np.abs(np.fft.rfft(y)) ** 2) / n          # periodogram
    f = np.fft.rfftfreq(n, d=1.0)
    P = P[1:]; f = f[1:]                            # drop the DC term
    hi = f > f_cut
    if hi.sum() < 4:
        return np.nan
    noise_psd = np.median(P[hi]) / LN2              # <-- exponential-median correction
    return float(np.clip(1.0 - noise_psd * len(P) / P.sum(), 0.0, 1.0))


# ------------------------------- validation ------------------------------- #
def _nanmedian(a):
    """np.nanmedian without the all-NaN RuntimeWarning (returns NaN for an empty/all-NaN
    slice, which is the intended sentinel)."""
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return float(np.median(a)) if a.size else np.nan


def validate(n_rep=200, seed=7, verbose=True):
    """Recover a KNOWN smooth signal fraction from synthetic (gap-free) data.

    Three signal shapes, four true SSF levels. Reports absolute bias per estimator.
    Reproduces the validation table in Methods 2.3 / Figure 5(a).
    """
    rng = np.random.default_rng(seed)
    t = np.arange(90)
    shapes = {
        "smooth sinusoid (28d)": np.sin(2 * np.pi * t / 28),
        "two-peak (like E2)":    np.exp(-((t % 28 - 13) ** 2) / 8)
                                 + 0.55 * np.exp(-((t % 28 - 21) ** 2) / 24),
        "sharp pulse (cramps)":  (np.minimum(t % 28, 28 - (t % 28)) < 2.5).astype(float),
    }
    err = {"ar1": [], "lin": [], "spec": []}
    if verbose:
        print(f"{'signal':>22} {'TRUE':>6} {'AR(1)':>8} {'ACF-lin':>9} {'SPECTRAL':>10}")
        print("-" * 60)
    for name, s in shapes.items():
        s = (s - s.mean()) / s.std()
        for true in (0.30, 0.40, 0.55, 0.70):
            est = {"ar1": [], "lin": [], "spec": []}
            sd = np.sqrt((1 - true) / true)
            for _ in range(n_rep):
                y = s + rng.normal(0, sd, len(s))
                est["ar1"].append(ssf_ar1(y))
                est["lin"].append(ssf_acf_linear(y))
                est["spec"].append(ssf_spectral(y))
            a, l, sp = (_nanmedian(est[k]) for k in ("ar1", "lin", "spec"))
            for k, v in zip(("ar1", "lin", "spec"), (a, l, sp)):
                err[k].append(abs(v - true))
            if verbose and true in (0.40, 0.70):
                print(f"{name:>22} {true:>6.2f} {a:>8.2f} {l:>9.2f} {sp:>10.2f}")
    if verbose:
        print("-" * 60)
        print(f"{'MEAN |BIAS|':>22} {'':>6} {np.mean(err['ar1']):>8.3f} "
              f"{np.mean(err['lin']):>9.3f} {np.mean(err['spec']):>10.3f}")
        print(f"{'MAX |BIAS|':>22} {'':>6} {np.max(err['ar1']):>8.3f} "
              f"{np.max(err['lin']):>9.3f} {np.max(err['spec']):>10.3f}")
    return err


def validate_gaps(n_rep=300, seed=9, missing=0.08, length=240, true=0.55, verbose=True):
    """Exercise the path the table above never touched: GAPPED series.

    Induce NaN gaps on a regular grid and compare the LEGACY 'compact' policy (which the
    pre-audit code used implicitly) against the corrected 'longest' policy. 'compact'
    is biased because it compresses the time axis across gaps; 'longest' is unbiased but
    honestly returns NaN when no contiguous run reaches the minimum length. Reporting the
    NaN-rate is part of the point: the correct method declines rather than fabricates.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(length)
    s = np.sin(2 * np.pi * t / 28); s = (s - s.mean()) / s.std()
    sd = np.sqrt((1 - true) / true)
    comp, longr = [], []
    for _ in range(n_rep):
        y = s + rng.normal(0, sd, length)
        y[rng.random(length) < missing] = np.nan            # NaN-holed regular grid
        comp.append(ssf_spectral(y, on_gap="compact"))
        longr.append(ssf_spectral(y, on_gap="longest"))
    comp = np.asarray(comp, float); longr = np.asarray(longr, float)
    out = dict(true=true, missing=missing,
               compact_median=_nanmedian(comp),
               longest_median=_nanmedian(longr),
               longest_nan_rate=float(np.mean(~np.isfinite(longr))),
               compact_bias=_nanmedian(comp) - true,
               longest_bias=_nanmedian(longr) - true)
    if verbose:
        print("\nGAP ROBUSTNESS (sinusoid, true SSF = %.2f, %d%% missing)" % (true, int(missing * 100)))
        print("-" * 60)
        print(f"  legacy 'compact' median : {out['compact_median']:.3f}   "
              f"(bias {out['compact_bias']:+.3f})")
        print(f"  corrected 'longest' med : {out['longest_median']:.3f}   "
              f"(bias {out['longest_bias']:+.3f}; NaN-rate {out['longest_nan_rate']:.2f})")
        print("  => 'compact' absorbs the gap into a biased number; 'longest' either")
        print("     estimates from a clean run or declines. It never silently compacts.")
    return out


if __name__ == "__main__":
    validate()
    validate_gaps()
