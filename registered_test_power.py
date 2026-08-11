#!/usr/bin/env python3
"""
Power of the EXACT preregistered test: phase-randomized surrogates.

The registered H3 test was NOT the LRT. It was a surrogate test for genuine
between-person heterogeneity in within-person coupling, with the predictor
phase-randomized (preserving the power spectrum, hence the full autocorrelation
structure, while destroying temporal alignment with the outcome).

Naive permutation is prohibited: it destroys autocorrelation, understates the null
variance of r, and inflates Type I error. (This correction was made before lodging.)

Test statistic:  S = SD of the within-person couplings r_i across individuals.
Null:            per-person phase randomization of the predictor.
p:               (1 + #{S_b >= S_obs}) / (B + 1)

RELIABILITY vs SSF -- READ THIS (audit note)
--------------------------------------------
The injected R_X / R_Y are RELIABILITIES: measurement error is scaled to the true-score
variance (signal + intrinsic white state fluctuation), exactly as in h4_frontier. They
are NOT smooth-signal fractions. This matters because:

  * In this generative model the white state noise (SIGMA_STATE) dominates the tiny
    cyclic signal (measured var_state / var_signal ~ 16), so the OUTCOME's actual
    smooth-signal fraction is ~0.03-0.07 -- an order of magnitude below the value one
    injects as R_Y. The predictor carries no state noise, so ITS SSF does match R_X.
  * The cycle coupling is attenuated by sqrt(SSF), not sqrt(reliability). Because the
    simulated outcome's SSF (~0.02 at R_Y=0.323) is far below the 0.323 the empirical
    section (7.3) treats as the outcome SSF, THIS SIMULATION MODELS A MUCH NOISIER
    OUTCOME THAN THE REAL INSTRUMENT, and therefore UNDERSTATES power. The reported
    power is a LOWER BOUND; the "SSF = 0.323" label on the manuscript's power table is
    a reliability, not an SSF.
  * Reconciling this is an AUTHORIAL decision, not a code patch: either (A) relabel the
    power table as reliability and declare the implied SSF, or (B) rebuild the generator
    to decouple coupling magnitude from smooth fraction so an outcome with SSF = 0.323
    can be simulated without inflating |r_i| past Table 1. A previous revision of this
    file tried to force the outcome SSF to 0.323; that target is UNREACHABLE here (state
    noise alone caps it near 0.066), so that change was reverted.

The diagnostic printed at the top of the power run reports the ACTUAL median SSF of the
simulated predictor and outcome, so the gap is visible rather than hidden.
"""
import argparse
import time
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

import h4_frontier as H
from ssf_estimators import ssf_spectral

# --- design and instruments as measured from mcPHASES ---
N_SUBJ    = 42
N_OBS     = 85      # median paired obs per person (corrected 6-level ordinal map)
SPAN      = 90      # days in study (~3 cycles)
# NB: injected as RELIABILITIES (see docstring). Set to the MEASURED SSF values, which is
# a CONSERVATIVE lower bound on reliability -> power reported here is a lower bound.
R_X_REAL  = 0.469   # predictor reliability (== measured E3G SSF; predictor has no state noise)
R_Y_REAL  = 0.323   # outcome reliability   (== measured self-report SSF; state-dominated)
R_X_IDEAL = 0.90    # "adequate instrument" reference for the ideal column (declared, not fitted)
R_Y_IDEAL = 0.90
B_SURR    = 500     # surrogates, as registered
N_REP     = 200     # replicates, as reported for Table 13
ALPHA     = 0.05


def phase_randomize(x, rng):
    """Preserve the power spectrum (hence autocorrelation); randomize phases."""
    n = len(x)
    X = np.fft.rfft(x)
    mag = np.abs(X)
    ph = rng.uniform(0, 2 * np.pi, len(X))
    ph[0] = 0.0                      # keep the DC term real
    if n % 2 == 0:
        ph[-1] = 0.0                 # keep Nyquist real
    return np.fft.irfft(mag * np.exp(1j * ph), n)


def _corr(a, b):
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def simulate_study(rng, sigma_b, R_x, R_y, n_subj=N_SUBJ, n_obs=N_OBS, span=SPAN):
    """Each person gets her OWN cycle phase offset and her own sampling days.

    Injection: classical RELIABILITY (measurement error scaled to true-score variance).
    State noise is part of the true score, as in h4_frontier. See module docstring for
    why the resulting outcome SSF is far below R_y.
    """
    b = rng.normal(H.DA_OPT - H.K_GAIN * H._E2_MEAN, sigma_b, n_subj)
    offs = rng.uniform(0, H.CYCLE_LEN, n_subj)          # individual cycle phase
    Xs, Ys = [], []
    for i in range(n_subj):
        days = np.sort(rng.choice(np.arange(span), size=n_obs, replace=False)).astype(float)
        x_true = H.e2(days + offs[i])                                  # smooth predictor, no state
        sig = H.inverted_u(b[i] + H.K_GAIN * x_true)                  # smooth outcome component
        y_true = sig + rng.normal(0, H.SIGMA_STATE, n_obs)            # + intrinsic state noise
        sdy = np.sqrt(max(y_true.var(), 1e-12) * (1 - R_y) / R_y)     # reliability injection
        sdx = np.sqrt(max(x_true.var(), 1e-12) * (1 - R_x) / R_x)
        Ys.append(y_true + rng.normal(0, sdy, n_obs))
        Xs.append(x_true + rng.normal(0, sdx, n_obs))
    return Xs, Ys, b


def surrogate_test(Xs, Ys, rng, B=B_SURR):
    """Registered test. Returns (p_value, S_obs, median S_null)."""
    r_obs = np.array([_corr(x, y) for x, y in zip(Xs, Ys)])
    S_obs = float(np.std(r_obs))
    S_null = np.empty(B)
    for k in range(B):
        r_s = np.array([_corr(phase_randomize(x, rng), y) for x, y in zip(Xs, Ys)])
        S_null[k] = np.std(r_s)
    p = (1 + int(np.sum(S_null >= S_obs))) / (B + 1)
    return p, S_obs, float(np.median(S_null))


def _ssf_diagnostic(rng):
    """Report the ACTUAL smooth-signal fraction produced by the injected reliabilities,
    so the reliability-vs-SSF gap is visible in the output."""
    Xs, Ys, _ = simulate_study(rng, 0.10, R_X_REAL, R_Y_REAL)
    ssx = float(np.nanmedian([ssf_spectral(np.asarray(x)) for x in Xs]))
    ssy = float(np.nanmedian([ssf_spectral(np.asarray(y)) for y in Ys]))
    print("DIAGNOSTIC (reliability vs SSF):")
    print(f"  injected reliabilities: Rx={R_X_REAL}, Ry={R_Y_REAL}")
    print(f"  ACTUAL median SSF of simulated series: predictor={ssx:.3f}, outcome={ssy:.3f}")
    print(f"  -> predictor SSF matches R_x (no state noise); outcome SSF << R_y because")
    print(f"     white state noise dominates. Cycle coupling attenuates by sqrt(SSF), so")
    print(f"     the power below is a LOWER BOUND vs a true-SSF-{R_Y_REAL} outcome.\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", default="power", choices=["power", "calib"])
    ap.add_argument("--reps", type=int, default=N_REP,
                    help="replicates per sigma_b (default 200, as reported for Table 13)")
    ap.add_argument("--surr", type=int, default=B_SURR,
                    help="phase surrogates per replicate (default 500, as registered)")
    ap.add_argument("--seed", type=int, default=23)
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    t0 = time.time()

    if a.mode == "calib":
        print("TYPE I CALIBRATION of the registered surrogate test")
        print("(sigma_b = 0: no coupling heterogeneity whatsoever)\n")
        n_rep = min(a.reps, 200)
        rej = 0
        for _ in range(n_rep):
            Xs, Ys, _ = simulate_study(rng, 0.0001, R_X_REAL, R_Y_REAL)
            p, _, _ = surrogate_test(Xs, Ys, rng, B=a.surr)
            rej += (p < ALPHA)
        fpr = rej / n_rep
        print(f"  false-positive rate = {fpr:.3f}   (nominal alpha = {ALPHA})")
        print("  -> the registered test is correctly calibrated." if abs(fpr - ALPHA) < 0.04
              else "  -> WARNING: miscalibrated; investigate before using.")

    else:
        print("POWER OF THE REGISTERED TEST (phase-randomized surrogates)")
        print(f"design: N={N_SUBJ}, {N_OBS} obs/person over {SPAN} days | "
              f"config: {a.reps} reps x {a.surr} surrogates\n")
        _ssf_diagnostic(rng)
        print(f"{'sigma_b':>8} {'power (actual instr.)':>22} {'power (ideal instr.)':>22}")
        print("-" * 56)
        rows = []
        grid = {"actual": (R_X_REAL, R_Y_REAL), "ideal": (R_X_IDEAL, R_Y_IDEAL)}
        for sb in [0.05, 0.075, 0.10, 0.15, 0.20, 0.30]:
            res = {}
            for tag, (rx, ry) in grid.items():
                rej = 0
                for _ in range(a.reps):
                    Xs, Ys, _ = simulate_study(rng, sb, rx, ry)
                    p, _, _ = surrogate_test(Xs, Ys, rng, B=a.surr)
                    rej += (p < ALPHA)
                res[tag] = rej / a.reps
            print(f"{sb:>8.3f} {res['actual']:>22.2f} {res['ideal']:>22.2f}", flush=True)
            rows.append(dict(sigma_b=sb, power_actual=res["actual"], power_ideal=res["ideal"]))
        pd.DataFrame(rows).to_csv("/mnt/user-data/outputs/registered_test_power.csv", index=False)
        print()
        print("  mcPHASES observed: fatigue surrogate p = 0.66 (null).")
        print("  The 'actual' column is a LOWER BOUND on power (see diagnostic): reliability")
        print("  was set to the measured SSF, and the outcome's true SSF is lower still.")
        print("  The sigma_b at which power crosses .80 is the largest heterogeneity the")
        print("  study could have MISSED -- an UPPER bound on that, given the lower-bound power.")

    print(f"\n[{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
