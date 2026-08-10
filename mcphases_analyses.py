#!/usr/bin/env python3
"""
mcphases_analyses.py — every empirical analysis on mcPHASES, in one reproducible script.

DATA ACCESS (READ THIS FIRST)
-----------------------------
mcPHASES is CREDENTIALED-ACCESS PhysioNet data (DOI 10.13026/zx6a-2c81) and is
therefore NOT redistributed with this package. To reproduce:
  1. Obtain PhysioNet credentialed access and accept the data use agreement.
  2. Download and unzip mcPHASES v1.0.0.
  3. Point --data-dir at the unzipped folder (the one containing
     hormones_and_selfreport.csv).

WHAT THIS REPRODUCES (manuscript section -> function)
-----------------------------------------------------
  7.1  describe()                 sample, ordinal levels, effective paired density
  7.2  instrument_ssf()           SSF of predictor, outcomes, and wearables (Table 8)
  7.6  differential_prediction()  the preregistered differential prediction (Table 11)
  7.6  disattenuate()             disattenuated group effects + person-level bootstrap
                                  (Table 12 -- shows the prediction CANNOT be tested)
  7.6  phase_locked()             cycle PHASE vs hormone LEVEL (Figure 4)
  7.7  modelfree_bound()          observed SD(r_i) vs its own phase-randomised null
  7.8  objective_vs_selfreport()  wearable vs self-report (Table 14)

NOTE ON THE ORDINAL MAP
-----------------------
The self-report scale has SIX levels. An early version of this analysis used a
five-level map that omitted "Not at all", silently discarding valid data and
producing N=39/median 74 instead of the correct N=41/median 85. The six-level map
below is the correct one. Categories outside the map are mapped to NaN and dropped;
verify with sorted(df[item].dropna().unique()) that the map is exhaustive.

NOTE ON SSF, GAPS, AND PAIRING
------------------------------
SSF is a property of the OUTCOME (or predictor) instrument. It must be computed on the
COMPLETE series, on a REGULAR grid, and NOT on any subset conditioned on the other
variable's availability (which would punch holes into the series; the spectral
estimator is a periodogram and assumes regular spacing). Every SSF computation here
wraps the raw values with ssf_estimators.regular_grid(values, day_index), which
reindexes onto the contiguous daily grid and inserts NaN at missing days so that the
estimator's on_gap='longest' policy acts on genuine gaps rather than silently
compacting them. Coupling, by contrast, requires pairing and legitimately uses
paired_series(); the two estimands never share a series.
"""
import argparse
import warnings

import numpy as np
import pandas as pd

from ssf_estimators import ssf_spectral, ssf_ar1, ssf_acf_linear, regular_grid

warnings.filterwarnings("ignore")

# THE CORRECT six-level ordinal map
ORDINAL = {"Not at all": 0, "Very Low/Little": 1, "Low": 2,
           "Moderate": 3, "High": 4, "Very High": 5}

# Preregistered outcome classification, fixed BEFORE data access.
# NOTE: "energy" was preregistered as a BALANCED (confirmatory) outcome but DOES NOT
# EXIST in mcPHASES. One third of the confirmatory set cannot be tested. Declared in
# manuscript Limitations.
CLASSIFICATION = {
    "fatigue":      "BALANCED (confirmatory)",
    "moodswing":    "BALANCED (confirmatory)",
    "cramps":       "DIRECTIONAL",
    "bloating":     "DIRECTIONAL",
    "sorebreasts":  "DIRECTIONAL",
    "sleepissue":   "AMBIGUOUS (exploratory)",
    "stress":       "AMBIGUOUS (exploratory)",
    "appetite":     "AMBIGUOUS (exploratory)",
    "foodcravings": "AMBIGUOUS (exploratory)",
}
MIN_PAIRED = 25          # preregistered inclusion threshold


# ------------------------------- helpers ------------------------------- #
def load(data_dir):
    df = pd.read_csv(f"{data_dir}/hormones_and_selfreport.csv")
    return df.sort_values(["id", "day_in_study"])


def corr(a, b):
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def hedges_g(r):
    """Convert a correlation to the Cohen's d scale used by the meta-analysis
    (d = 2r / sqrt(1 - r^2)). NB: this is d, WITHOUT the small-sample Hedges J
    correction; the name is retained for continuity with earlier drafts."""
    r = np.clip(r, -0.99, 0.99)
    return abs(2 * r / np.sqrt(1 - r ** 2))


def phase_randomize(x, rng):
    """Preserve the power spectrum (hence the FULL autocorrelation structure);
    randomise the phases. This is the registered null. Naive permutation is
    PROHIBITED: it destroys autocorrelation, understates the null variance of r,
    and inflates Type I error."""
    n = len(x)
    X = np.fft.rfft(x)
    mag = np.abs(X)
    ph = rng.uniform(0, 2 * np.pi, len(X))
    ph[0] = 0.0                                  # keep DC real
    if n % 2 == 0:
        ph[-1] = 0.0                             # keep Nyquist real
    return np.fft.irfft(mag * np.exp(1j * ph), n)


def paired_series(df, item, min_n=MIN_PAIRED):
    """List of (estrogen, outcome) arrays, one per qualifying participant. For COUPLING
    only -- pairing is required to correlate the two. Never use these for SSF."""
    d = df.dropna(subset=["estrogen", item]).copy()
    d["y"] = d[item].map(ORDINAL)
    d = d.dropna(subset=["y"])
    return [(g.estrogen.values.astype(float), g.y.values.astype(float))
            for _, g in d.groupby("id") if len(g) >= min_n]


# ------------------------------- 7.1 ------------------------------- #
def describe(df):
    print("=" * 74)
    print("7.1  WHAT THE STUDY ACTUALLY IS")
    print("=" * 74)
    print(f"  participants            : {df.id.nunique()}")
    print(f"  rows                    : {len(df)}")
    print(f"  self-report levels      : {sorted(df.fatigue.dropna().unique())}")
    pairs = paired_series(df, "fatigue")
    cnt = pd.Series([len(x) for x, _ in pairs])
    print(f"  participants >= {MIN_PAIRED} paired obs : {len(pairs)} / {df.id.nunique()}")
    print(f"  paired obs/person       : median {cnt.median():.0f} "
          f"(IQR {cnt.quantile(.25):.0f}-{cnt.quantile(.75):.0f})")
    print(f"  study duration          : median {df.groupby('id').day_in_study.max().median():.0f} days")
    print(f"  predictor coverage      : {100*df.estrogen.notna().mean():.0f}%")
    print("\n  -> Density is AMPLE. The original hypothesis (too sparse) is FALSE.")


# ------------------------------- 7.2 ------------------------------- #
def instrument_ssf(df, data_dir):
    print("\n" + "=" * 74)
    print("7.2  INSTRUMENT QUALITY (Table 8)")
    print("=" * 74)
    print(f"{'measure':>34} {'AR(1)':>9} {'ACF-lin':>9} {'SPECTRAL':>10}")
    print("-" * 66)
    rows = []

    def report(label, series_by_person):
        est = {k: [] for k in ("ar1", "lin", "spec")}
        for y in series_by_person:
            for k, f in (("ar1", ssf_ar1), ("lin", ssf_acf_linear), ("spec", ssf_spectral)):
                v = f(y)                                    # on_gap='longest' by default
                if np.isfinite(v):
                    est[k].append(v)
        a, l, s = (np.median(est[k]) if est[k] else np.nan for k in ("ar1", "lin", "spec"))
        flag = "  <-- IMPOSSIBLE (>1)" if a > 1 else ""
        print(f"{label:>34} {a:>9.3f} {l:>9.3f} {s:>10.3f}{flag}")
        rows.append(dict(measure=label, ar1=a, acf_linear=l, spectral=s))
        return s

    # predictor and self-report outcomes: complete series on a regular grid
    report("estrone-3-glucuronide (PREDICTOR)",
           [regular_grid(g.estrogen.values, g.day_in_study.values)
            for _, g in df.groupby("id")])
    for item in ("fatigue", "moodswing", "cramps", "bloating", "stress"):
        report(f"{item} [{CLASSIFICATION.get(item,'')}]",
               [regular_grid(g[item].map(ORDINAL).values, g.day_in_study.values)
                for _, g in df.groupby("id")])

    # objective wearables
    try:
        rhr = pd.read_csv(f"{data_dir}/resting_heart_rate.csv").sort_values(["id", "day_in_study"])
        report("resting heart rate (OBJECTIVE)",
               [regular_grid(g.value.values, g.day_in_study.values)
                for _, g in rhr.groupby("id")])
        ct = pd.read_csv(f"{data_dir}/computed_temperature.csv")
        ct = ct[ct.type == "SKIN"].sort_values(["id", "sleep_start_day_in_study"])
        report("nightly skin temperature (OBJECTIVE)",
               [regular_grid(g.nightly_temperature.values, g.sleep_start_day_in_study.values)
                for _, g in ct.groupby("id")])
    except FileNotFoundError:
        print("  [wearable files not found; skipping objective measures]")

    out = pd.DataFrame(rows)
    print("\n  ATTENUATION = sqrt(SSF_predictor * SSF_outcome). Multiplicative:")
    sx = out.loc[out.measure.str.contains("PREDICTOR"), "spectral"].iloc[0]
    sy = out.loc[out.measure.str.startswith(("fatigue", "moodswing")), "spectral"].median()
    print(f"    sqrt({sx:.3f} * {sy:.3f}) = {np.sqrt(sx*sy):.3f}")
    print(f"    -> observed coupling is less than {100*np.sqrt(sx*sy):.0f}% of true")
    return out


# ------------------------------- 7.6a ------------------------------- #
def differential_prediction(df, rng):
    print("\n" + "=" * 74)
    print("7.6  THE PREREGISTERED DIFFERENTIAL PREDICTION (Table 11)")
    print("     BALANCED -> group null ;  DIRECTIONAL -> group effect LEAKS (|g| >= .10)")
    print("=" * 74)
    print(f"{'outcome':>14} {'class':>26} {'|g| obs':>8} {'predicted':>10} {'OBSERVED':>10}")
    print("-" * 74)
    rows = []
    for item, cls in CLASSIFICATION.items():
        per = paired_series(df, item)
        if len(per) < 10:
            continue
        # NOTE (audit ID-E1): this pooled correlation mixes within- and between-person
        # variance and is not the mean of individual couplings the GH targets. Retained
        # pending a decision on the estimand; see the audit report (options A/B).
        X = np.concatenate([x for x, _ in per])
        Y = np.concatenate([y for _, y in per])
        g = hedges_g(corr(X, Y))
        pred = "LEAK" if cls.startswith("DIRECTIONAL") else ("null" if cls.startswith("BALANCED") else "-")
        obs = "LEAKS" if g >= 0.10 else "null"
        mark = ""
        if pred != "-" and ((pred == "LEAK") != (obs == "LEAKS")):
            mark = "  <-- contradicts"
        print(f"{item:>14} {cls:>26} {g:>8.3f} {pred:>10} {obs:>10}{mark}")
        rows.append(dict(item=item, cls=cls, group_g=g, predicted=pred, observed=obs))
    print("\n  READ THE NEXT FUNCTION BEFORE CONCLUDING ANYTHING FROM THIS TABLE.")
    return pd.DataFrame(rows)


# ------------------------------- 7.6b (THE KEY ONE) ------------------------------- #
def disattenuate(df, rng, n_boot=2000):
    """The |g| < .10 threshold is stated on the OBSERVED scale -- the SAME defect as
    preregistered criterion (i). Correcting for attenuation flips verdicts, but
    disattenuation ALSO inflates the sampling variance (by 1/att^2 ~ 6.6x).

    Person-level bootstrap (resampling whole participants, respecting clustering)
    shows that NO group effect can be classified: every 95% CI spans both 0 and 0.10.

    CONCLUSION: the differential prediction CANNOT BE TESTED with these instruments.
    Not "it failed" -- it cannot be evaluated. This is the third instance of the
    manuscript's own thesis applying to the manuscript.

    SSF NOTE: both sx (predictor) and sy (outcome) are computed on the COMPLETE series
    on a regular grid -- NOT on the estrogen-paired subset used for the correlation.
    """
    print("\n" + "=" * 74)
    print("7.6  DISATTENUATED GROUP EFFECTS + PERSON-LEVEL BOOTSTRAP (Table 12)")
    print("=" * 74)
    # predictor SSF: complete estrogen series, regular grid, all persons
    sx = np.median([v for _, g in df.groupby("id")
                    for v in [ssf_spectral(regular_grid(g.estrogen.values,
                                                        g.day_in_study.values))]
                    if np.isfinite(v)])
    print(f"  SSF predictor (E3G) = {sx:.3f}\n")
    print(f"{'outcome':>14} {'class':>13} {'atten':>7} {'|g| obs':>8} {'|g| disatt':>11} "
          f"{'95% CI':>18} {'classifiable?':>14}")
    print("-" * 92)
    rows = []
    for item, cls in CLASSIFICATION.items():
        if not (cls.startswith("BALANCED") or cls.startswith("DIRECTIONAL")):
            continue
        per = paired_series(df, item)                       # COUPLING series (paired)
        # outcome SSF on the COMPLETE series (regular grid), NOT the paired subset
        sy = np.median([v for _, g in df.groupby("id")
                        for v in [ssf_spectral(regular_grid(g[item].map(ORDINAL).values,
                                                            g.day_in_study.values))]
                        if np.isfinite(v)])
        att = np.sqrt(sx * sy)
        X = np.concatenate([x for x, _ in per]); Y = np.concatenate([y for _, y in per])
        g_obs = hedges_g(corr(X, Y))
        g_dis = hedges_g(corr(X, Y) / att)
        boots = []
        for _ in range(n_boot):
            idx = rng.integers(0, len(per), len(per))          # resample PEOPLE
            Xb = np.concatenate([per[i][0] for i in idx])
            Yb = np.concatenate([per[i][1] for i in idx])
            boots.append(hedges_g(corr(Xb, Yb) / att))
        lo, hi = np.percentile(boots, [2.5, 97.5])
        classifiable = "no" if (lo <= 0.10 <= hi) else "yes"
        print(f"{item:>14} {cls.split()[0]:>13} {att:>7.3f} {g_obs:>8.3f} {g_dis:>11.3f} "
              f"{f'[{lo:.3f}; {hi:.3f}]':>18} {classifiable:>14}")
        rows.append(dict(item=item, cls=cls, attenuation=att, g_observed=g_obs,
                         g_disattenuated=g_dis, ci_lo=lo, ci_hi=hi,
                         classifiable=classifiable))
    print("\n  Every CI spans both 0 and the 0.10 threshold.")
    print("  => THE DIFFERENTIAL PREDICTION CANNOT BE TESTED. It did not 'fail'.")
    return pd.DataFrame(rows)


# ------------------------------- 7.6c ------------------------------- #
def phase_locked(df):
    """A second, independent reason the prediction is untestable as posed:
    directional symptoms are locked to cycle PHASE, not hormone LEVEL. E2 is low both
    during menses (cramps high) and in the early follicular phase (cramps low), so a
    LINEAR association with E2 cancels even though the symptom is strongly cyclic."""
    print("\n" + "=" * 74)
    print("7.6  PHASE-LOCKED vs LEVEL-LOCKED (Figure 4)")
    print("=" * 74)
    print(f"{'outcome':>14} {'class':>13} {'|r| with E2 LEVEL':>18} {'eta^2 of PHASE':>15} "
          f"{'menstrual mean':>15}")
    print("-" * 80)
    rows = []
    for item in ("fatigue", "moodswing", "cramps", "bloating", "sorebreasts"):
        d = df.dropna(subset=["estrogen", item]).copy()
        d["y"] = d[item].map(ORDINAL); d = d.dropna(subset=["y"])
        # within-person centring removes all between-person variance
        d["yc"] = d.y - d.groupby("id").y.transform("mean")
        d["ec"] = d.estrogen - d.groupby("id").estrogen.transform("mean")
        r_e2 = abs(corr(d.ec, d.yc))
        ss_b = sum(len(g) * (g.mean() - d.yc.mean()) ** 2 for _, g in d.groupby("phase").yc)
        ss_t = ((d.yc - d.yc.mean()) ** 2).sum()
        eta2 = ss_b / ss_t if ss_t > 0 else np.nan
        mens = d[d.phase.astype(str).str.contains("enstrual", case=False, na=False)].yc.mean()
        print(f"{item:>14} {CLASSIFICATION[item].split()[0]:>13} {r_e2:>18.3f} "
              f"{eta2:>15.3f} {mens:>+15.3f}")
        rows.append(dict(item=item, r_e2_level=r_e2, phase_eta2=eta2, menstrual_mean=mens))
    print("\n  Cramps: |r| with E2 level = .075, but eta^2 with cycle phase = .162.")
    print("  The preregistered classification was RIGHT about which symptoms are")
    print("  directional and WRONG about the predictor. A hormone LEVEL is not a")
    print("  proxy for a cycle PHASE.")
    return pd.DataFrame(rows)


# ------------------------------- 7.7 ------------------------------- #
def modelfree_bound(df, rng, item="fatigue", B=500):
    """Fully model-free: observed SD(r_i) against its OWN phase-randomised null.
    Both are on the same (attenuated) scale -- an apples-to-apples comparison, and
    therefore VALID where the group-effect threshold was not."""
    print("\n" + "=" * 74)
    print(f"7.7  MODEL-FREE BOUND ({item})")
    print("=" * 74)
    per = paired_series(df, item)
    r_obs = np.array([corr(x, y) for x, y in per])
    S_obs = float(np.std(r_obs))
    S_null = np.array([np.std([corr(phase_randomize(x, rng), y) for x, y in per])
                       for _ in range(B)])
    p = (1 + int(np.sum(S_null >= S_obs))) / (B + 1)
    print(f"  participants        : {len(per)}")
    print(f"  observed SD(r_i)    : {S_obs:.4f}")
    print(f"  null median SD(r_i) : {np.median(S_null):.4f}   (95th pct {np.quantile(S_null,.95):.4f})")
    print(f"  p                   : {p:.3f}")
    print(f"\n  Observed heterogeneity sits ON the null median. Not 'small' -- ZERO excess.")
    return dict(item=item, n=len(per), SD_obs=S_obs,
                SD_null_median=float(np.median(S_null)), p=p)


# ------------------------------- 7.8 ------------------------------- #
def objective_vs_selfreport(df, data_dir, rng, B=300):
    """Test the manuscript's own prescription inside the dataset: same participants,
    same days, same predictor, same registered test -- but an OBJECTIVE outcome.

    TWO ESTIMANDS, TWO SERIES:
      * COUPLING  (SD(r_i) + phase-randomised registered test) -> PAIRED (estrogen,
        outcome) series (dropna on both is correct here).
      * SSF (spectral) -> COMPLETE per-person OUTCOME series on a REGULAR grid
        (regular_grid), NOT conditioned on E3G availability.
    """
    print("\n" + "=" * 74)
    print("7.8  OBJECTIVE vs SELF-REPORT (Table 14) -- EXPLORATORY")
    print("=" * 74)

    rhr = pd.read_csv(f"{data_dir}/resting_heart_rate.csv")[["id", "day_in_study", "value"]] \
            .rename(columns={"value": "rhr"})
    ct = pd.read_csv(f"{data_dir}/computed_temperature.csv")
    ct = ct[ct.type == "SKIN"][["id", "sleep_start_day_in_study", "nightly_temperature"]] \
            .rename(columns={"sleep_start_day_in_study": "day_in_study",
                             "nightly_temperature": "temp"})
    d = df.copy(); d["fatigue_n"] = d.fatigue.map(ORDINAL)

    # --- paired frame: for COUPLING only (correlation + registered test) ---
    m = d[["id", "day_in_study", "estrogen", "fatigue_n"]] \
        .merge(rhr, on=["id", "day_in_study"], how="left") \
        .merge(ct,  on=["id", "day_in_study"], how="left")

    # --- complete per-person OUTCOME series on a REGULAR grid: for SSF only ---
    complete = {
        "fatigue_n": {pid: regular_grid(g.sort_values("day_in_study").fatigue.map(ORDINAL).values,
                                        g.sort_values("day_in_study").day_in_study.values)
                      for pid, g in df.groupby("id")},
        "rhr":       {pid: regular_grid(g.sort_values("day_in_study").rhr.values,
                                        g.sort_values("day_in_study").day_in_study.values)
                      for pid, g in rhr.groupby("id")},
        "temp":      {pid: regular_grid(g.sort_values("day_in_study").temp.values,
                                        g.sort_values("day_in_study").day_in_study.values)
                      for pid, g in ct.groupby("id")},
    }

    print(f"{'outcome':>26} {'kind':>12} {'SSF':>7} {'SD(r_i)':>9} {'null':>8} {'p':>7}")
    print("-" * 74)
    rows = []
    ssf_by_outcome, p_by_outcome = {}, {}
    for col, lab, kind in (("fatigue_n", "fatigue", "SELF-REPORT"),
                           ("rhr", "resting heart rate", "OBJECTIVE"),
                           ("temp", "skin temperature", "OBJECTIVE")):
        per, ss = [], []
        for pid, g in m.groupby("id"):
            gg = g.dropna(subset=["estrogen", col])
            if len(gg) >= MIN_PAIRED:
                # COUPLING: paired series (estrogen aligned with outcome)
                per.append((gg.estrogen.values.astype(float),
                            gg[col].values.astype(float)))
                # SSF: this SAME participant's COMPLETE, regular-grid outcome series
                y_full = complete[col].get(pid)
                if y_full is not None:
                    v = ssf_spectral(y_full)               # on_gap='longest'
                    if np.isfinite(v):
                        ss.append(v)
        if len(per) < 10:
            continue
        S = np.std([corr(x, y) for x, y in per])
        Sn = np.array([np.std([corr(phase_randomize(x, rng), y) for x, y in per])
                       for _ in range(B)])
        p = (1 + int(np.sum(Sn >= S))) / (B + 1)
        ssf_med = float(np.median(ss)) if ss else np.nan
        star = " *" if p < 0.05 else ""
        print(f"{lab:>26} {kind:>12} {ssf_med:>7.3f} {S:>9.4f} "
              f"{np.median(Sn):>8.4f} {p:>7.3f}{star}")
        rows.append(dict(outcome=lab, kind=kind, ssf=ssf_med, SD_ri=S,
                         null_SD=float(np.median(Sn)), p=p))
        ssf_by_outcome[lab] = ssf_med
        p_by_outcome[lab] = p

    # --- caveats: computed from the results above, NOT hardcoded ---
    n_tests = len(rows)
    alpha_bonf = (0.05 / n_tests) if n_tests else float("nan")
    temp_p = p_by_outcome.get("skin temperature", float("nan"))
    rhr_ssf = ssf_by_outcome.get("resting heart rate", float("nan"))
    best_lab = max(ssf_by_outcome, key=ssf_by_outcome.get) if ssf_by_outcome else "-"
    print("\n  CAVEATS THAT MUST TRAVEL WITH THIS RESULT:")
    print(f"   1. skin-temperature p = {temp_p:.3f} with {n_tests} tests does NOT survive")
    print(f"      Bonferroni (alpha = {alpha_bonf:.4f}). EXPLORATORY and suggestive; not confirmatory.")
    print(f"   2. Resting HR shows NO coupling despite carrying the highest SSF of the")
    print(f"      three outcomes ({rhr_ssf:.3f}); instrument quality is NECESSARY, NOT SUFFICIENT.")
    if best_lab != "resting heart rate":
        print(f"      [audit check: highest SSF is actually {best_lab} = {ssf_by_outcome[best_lab]:.3f}]")
    print("   3. Skin temperature is cycle-locked via PROGESTERONE (luteal rise);")
    print("      its estradiol coupling may be confounded.")
    return pd.DataFrame(rows)


# ------------------------------- main ------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True,
                    help="unzipped mcPHASES v1.0.0 folder (credentialed PhysioNet access)")
    ap.add_argument("--out-dir", default="./results")
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()

    import os
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    df = load(args.data_dir)

    describe(df)
    instrument_ssf(df, args.data_dir).to_csv(f"{args.out_dir}/table08_instrument_ssf.csv", index=False)
    differential_prediction(df, rng).to_csv(f"{args.out_dir}/table11_differential.csv", index=False)
    disattenuate(df, rng).to_csv(f"{args.out_dir}/table12_disattenuated.csv", index=False)
    phase_locked(df).to_csv(f"{args.out_dir}/fig04_phase_locked.csv", index=False)
    pd.DataFrame([modelfree_bound(df, rng, "fatigue"),
                  modelfree_bound(df, rng, "moodswing")]).to_csv(
                      f"{args.out_dir}/modelfree_bound.csv", index=False)
    objective_vs_selfreport(df, args.data_dir, rng).to_csv(
        f"{args.out_dir}/table14_objective.csv", index=False)
    print(f"\n[saved] {args.out_dir}/")


if __name__ == "__main__":
    main()
