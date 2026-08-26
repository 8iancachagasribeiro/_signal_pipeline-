#!/usr/bin/env python3
"""
period_search_test.py — Does estimating the period per person rescue harmonic regression?

WHY THIS TEST COMES FIRST
-------------------------
Table 5 of the proposal shows harmonic regression degrading from |r| = 0.999 to 0.868 as
cycle-length variability rises to a standard deviation of 5 days, crossing below PCA
(flat at 0.887). That degradation is the empirical justification for the wavelet route,
because wavelets handle non-stationary period.

But there is a much cheaper fix for the same problem: stop fixing T = 28 and ESTIMATE the
period per person, by fitting the harmonic model across a grid of candidate periods and
keeping the best. This is what the Lomb-Scargle periodogram does natively. It costs a few
lines and no new method.

A reviewer who knows the field will ask exactly this: why wavelets, when you could just
estimate the period?

The answer must be measured, not argued:

  - If the period search recovers the loss, the wavelet justification WEAKENS and the
    proposal has to say so.
  - If it does not recover the loss, the wavelet justification becomes MEASURED rather
    than asserted, which is stronger than what the proposal currently claims.

Either outcome is publishable and either outcome is useful. Only running it settles it.

USAGE
-----
    python period_search_test.py
"""
import warnings

import numpy as np

warnings.filterwarnings("ignore")

T_DAYS = 90
N_SUB = 60
SEED = 2026
NOMINAL = 28.0

BASE = np.array([62.0, 34.5, 33.8, 40.0, 800.0, 600.0])
AMPL = np.array([2.5, 0.35, 0.30, 4.0, 60.0, 45.0])
SIGNS = np.array([+1, +1, +1, -1, -1, +1])
SSF_TARGET = np.array([0.47, 0.34, 0.34, 0.40, 0.35, 0.35])

# candidate periods: the physiological range of the menstrual cycle
PERIOD_GRID = np.arange(21.0, 35.5, 0.5)


def make_subject(rng, jitter):
    """Six non-negative channels sharing one latent cycle of VARIABLE duration."""
    t = np.arange(T_DAYS, dtype=float)
    edges = [0.0]
    while edges[-1] < T_DAYS:
        edges.append(edges[-1] + rng.normal(NOMINAL, jitter))
    phase = np.interp(t, edges, np.arange(len(edges)) * 2 * np.pi)
    lat = np.sin(phase)
    lat = (lat - lat.mean()) / lat.std()

    X = np.empty((T_DAYS, 6))
    for j in range(6):
        sig = SIGNS[j] * AMPL[j] * lat
        s = SSF_TARGET[j]
        sd = np.sqrt(max(sig.var(), 1e-12) * (1 - s) / s)
        X[:, j] = BASE[j] + sig + rng.normal(0, sd, T_DAYS)
    return np.clip(X, 1e-9, None), lat, t


def pca_recover(X):
    sd = X.std(0)
    Z = (X - X.mean(0)) / np.where(sd < 1e-9, 1.0, sd)
    Zc = Z - Z.mean(0)
    U, S, Vt = np.linalg.svd(Zc, full_matrices=False)
    return Zc @ Vt.T[:, 0]


def _design(t, period, n_harm):
    cols = [np.ones_like(t)]
    for h in range(1, n_harm + 1):
        cols += [np.cos(2 * np.pi * h * t / period),
                 np.sin(2 * np.pi * h * t / period)]
    return np.column_stack(cols)


def _fit_channel(y, t, period, n_harm):
    """Returns (phase of first harmonic, R^2 of the fit)."""
    D = _design(t, period, n_harm)
    beta, *_ = np.linalg.lstsq(D, y, rcond=None)
    resid = y - D @ beta
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1.0 - (resid ** 2).sum() / ss_tot if ss_tot > 0 else 0.0
    return np.arctan2(beta[2], beta[1]), float(np.clip(r2, 0.0, 0.999))


def harmonic_recover(X, t, n_harm=2, period=NOMINAL):
    """Harmonic regression at a GIVEN period."""
    phases, weights = [], []
    for j in range(X.shape[1]):
        ph, r2 = _fit_channel(X[:, j], t, period, n_harm)
        phases.append(ph)
        weights.append(r2 / (1.0 - r2))
    phases, weights = np.array(phases), np.array(weights)
    if weights.sum() <= 0:
        return np.zeros_like(t), 0.0
    ref = phases[np.argmax(weights)]
    aligned = phases + np.pi * (np.abs(np.angle(np.exp(1j * (phases - ref)))) > np.pi / 2)
    C = np.sum(weights * np.cos(aligned))
    S = np.sum(weights * np.sin(aligned))
    phi = np.arctan2(S, C)
    return np.cos(2 * np.pi * t / period - phi), float(weights.sum())


def harmonic_search(X, t, n_harm=2, grid=PERIOD_GRID):
    """Harmonic regression with the period ESTIMATED per person.

    The selection criterion is the total signal-to-noise of the multi-channel fit,
    summed across channels. This is the Lomb-Scargle logic: scan candidate periods,
    keep the one the data support best. No knowledge of the true period is used.
    """
    best_score, best_rec, best_p = -np.inf, None, None
    for p in grid:
        rec, score = harmonic_recover(X, t, n_harm, period=p)
        if score > best_score:
            best_score, best_rec, best_p = score, rec, p
    return best_rec, best_p


def main():
    print("=" * 70)
    print("DOES ESTIMATING THE PERIOD RESCUE HARMONIC REGRESSION?")
    print("=" * 70)
    print(f"  {N_SUB} subjects, {T_DAYS} days, nominal cycle {NOMINAL:.0f} days")
    print(f"  period grid searched: {PERIOD_GRID[0]:.0f} to {PERIOD_GRID[-1]:.0f} days,"
          f" step {PERIOD_GRID[1]-PERIOD_GRID[0]:.1f}")
    print()
    print(f"{'SD of cycle length':>19} {'PCA':>8} {'harmonic T=28':>15} "
          f"{'harmonic T searched':>21} {'median T found':>16}")
    print("-" * 84)

    for jit in (0.0, 1.5, 3.0, 5.0):
        rng = np.random.default_rng(SEED)
        r_pca, r_fix, r_srch, p_found = [], [], [], []
        for _ in range(N_SUB):
            X, lat, t = make_subject(rng, jit)
            r_pca.append(abs(np.corrcoef(pca_recover(X), lat)[0, 1]))
            rec_f, _ = harmonic_recover(X, t, 2, NOMINAL)
            r_fix.append(abs(np.corrcoef(rec_f, lat)[0, 1]))
            rec_s, p = harmonic_search(X, t, 2)
            r_srch.append(abs(np.corrcoef(rec_s, lat)[0, 1]))
            p_found.append(p)
        print(f"{jit:>17.1f}d {np.median(r_pca):>8.3f} {np.median(r_fix):>15.3f} "
              f"{np.median(r_srch):>21.3f} {np.median(p_found):>15.1f}d")

    print()
    print("=" * 70)
    print("HOW TO READ")
    print("=" * 70)
    print("  If 'T searched' recovers close to the 0.999 obtained at zero variability,")
    print("  the period search solves the problem and the WAVELET JUSTIFICATION WEAKENS.")
    print()
    print("  If it does not recover, the wavelet justification becomes MEASURED rather")
    print("  than asserted: the period is not merely unknown, it is non-stationary WITHIN")
    print("  a series, which a single global period cannot represent however well chosen.")


if __name__ == "__main__":
    main()
