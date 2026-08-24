#!/usr/bin/env python3
"""Audited mcPHASES empirical analyses for Article 1.

Audit principles
----------------
1. A participant may contribute more than one study_interval. Each interval is treated
   as a separate temporal segment for FFT/surrogate operations; intervals are never
   bridged across the long calendar gap.
2. Participant-level coupling is obtained by Fisher-z aggregation of qualifying
   interval-specific correlations, so repeated intervals do not become independent
   participants.
3. Pooled within-person associations are centered within participant-interval.
4. Fitbit files described as daily are collapsed to one value per participant-interval-day
   before merging. Resting-heart-rate values <= 0 are treated as missing because 0 bpm is
   physiologically impossible and the file README defines the field as daily resting bpm.
5. SSF is a smooth-signal fraction, not classical reliability. Any attenuation calculation
   using SSF is presented only as an identification sensitivity analysis.
6. Phase-randomized surrogates operate on the regular daily calendar within each interval,
   filling internal missing predictor days only to construct the FFT surrogate and then
   sampling the surrogate back at the originally observed days.
"""
from __future__ import annotations
import argparse, os, warnings
import numpy as np
import pandas as pd
from ssf_estimators import ssf_spectral, ssf_ar1, ssf_acf_linear, regular_grid

warnings.filterwarnings("ignore")

ORDINAL = {
    "Not at all": 0, "Very Low/Little": 1, "Low": 2,
    "Moderate": 3, "High": 4, "Very High": 5,
}
CLASSIFICATION = {
    "fatigue": "BALANCED (confirmatory)",
    "moodswing": "BALANCED (confirmatory)",
    "cramps": "DIRECTIONAL",
    "bloating": "DIRECTIONAL",
    "sorebreasts": "DIRECTIONAL",
    "sleepissue": "AMBIGUOUS (exploratory)",
    "stress": "AMBIGUOUS (exploratory)",
    "appetite": "AMBIGUOUS (exploratory)",
    "foodcravings": "AMBIGUOUS (exploratory)",
}
MIN_PAIRED = 25


def corr0(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b); a, b = a[ok], b[ok]
    if len(a) < 4 or np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return 0.0
    ac, bc = a - a.mean(), b - b.mean()
    return float(np.dot(ac, bc) / np.sqrt(np.dot(ac, ac) * np.dot(bc, bc)))


def effect_g(r):
    r = float(np.clip(r, -0.99, 0.99))
    return abs(2 * r / np.sqrt(1 - r * r))


def fisher_aggregate(rs, ns=None):
    rs = np.asarray(rs, float)
    if ns is None:
        ns = np.full(len(rs), 4.0)
    ns = np.asarray(ns, float)
    ok = np.isfinite(rs) & (ns > 3)
    if not ok.any():
        return np.nan
    z = np.arctanh(np.clip(rs[ok], -0.999999, 0.999999))
    w = np.maximum(ns[ok] - 3, 1.0)
    return float(np.tanh(np.sum(w * z) / np.sum(w)))


def phase_randomize(x, rng):
    x = np.asarray(x, float); n = len(x)
    X = np.fft.rfft(x); mag = np.abs(X)
    ph = rng.uniform(0, 2 * np.pi, len(X)); ph[0] = 0.0
    if n % 2 == 0:
        ph[-1] = 0.0
    return np.fft.irfft(mag * np.exp(1j * ph), n)


def prepare_calendar(x, day):
    """Return interpolated regular-grid predictor plus positions of observed days."""
    day = np.asarray(day, int); x = np.asarray(x, float)
    order = np.argsort(day, kind="mergesort"); day, x = day[order], x[order]
    if np.any(np.diff(day) == 0):
        d = pd.DataFrame({"day": day, "x": x}).groupby("day", as_index=False).x.mean()
        day, x = d.day.to_numpy(int), d.x.to_numpy(float)
    lo, hi = int(day.min()), int(day.max())
    pos = day - lo
    grid = np.full(hi - lo + 1, np.nan); grid[pos] = x
    finite = np.isfinite(grid)
    if not finite.all():
        z = np.arange(len(grid))
        grid[~finite] = np.interp(z[~finite], z[finite], grid[finite])
    return grid, pos


def load_hormones(data_dir):
    d = pd.read_csv(os.path.join(data_dir, "hormones_and_selfreport.csv"))
    d["episode"] = d["id"].astype(str) + "_" + d["study_interval"].astype(str)
    return d.sort_values(["id", "study_interval", "day_in_study"])


def item_records(df, item, min_n=MIN_PAIRED):
    d = df[["id", "study_interval", "episode", "day_in_study", "estrogen", item]].copy()
    d["y"] = d[item].map(ORDINAL)
    d = d.dropna(subset=["estrogen", "y", "day_in_study"])
    out = []
    for (pid, interval, episode), g in d.groupby(["id", "study_interval", "episode"]):
        g = g.groupby("day_in_study", as_index=False).agg(x=("estrogen", "mean"), y=("y", "mean"))
        g = g.sort_values("day_in_study")
        if len(g) < min_n:
            continue
        day = g.day_in_study.to_numpy(int)
        x = g.x.to_numpy(float); y = g.y.to_numpy(float)
        grid, pos = prepare_calendar(x, day)
        out.append(dict(id=int(pid), interval=interval, episode=episode,
                        x=x, y=y, day=day, n=len(g), grid=grid, pos=pos))
    return out


def participant_r(records, rng=None):
    by = {}
    for r in records:
        x = phase_randomize(r["grid"], rng)[r["pos"]] if rng is not None else r["x"]
        rv = corr0(x, r["y"])
        by.setdefault(r["id"], []).append((rv, r["n"]))
    return {
        pid: fisher_aggregate([r for r, n in vals], [n for r, n in vals])
        for pid, vals in by.items()
    }


def pooled_episode_centered(records):
    xs, ys = [], []
    for r in records:
        xs.append(r["x"] - np.mean(r["x"]))
        ys.append(r["y"] - np.mean(r["y"]))
    return corr0(np.concatenate(xs), np.concatenate(ys)) if xs else np.nan


def daily_rhr(data_dir):
    d = pd.read_csv(os.path.join(data_dir, "resting_heart_rate.csv"))
    d.loc[d["value"] <= 0, "value"] = np.nan
    return (d.groupby(["id", "study_interval", "day_in_study"], as_index=False)
              .value.median().dropna(subset=["value"]))


def daily_temp(data_dir):
    d = pd.read_csv(os.path.join(data_dir, "computed_temperature.csv"))
    d = d[d["type"] == "SKIN"]
    d = (d.groupby(["id", "study_interval", "sleep_start_day_in_study"], as_index=False)
           .nightly_temperature.mean())
    return d.rename(columns={"sleep_start_day_in_study": "day_in_study",
                             "nightly_temperature": "value"})


def participant_balanced_ssf(d, value_col, day_col="day_in_study"):
    rows = []
    for (pid, interval), g in d.groupby(["id", "study_interval"]):
        y = regular_grid(g[value_col].values, g[day_col].values)
        rows.append(dict(id=pid, interval=interval,
                         ar1=ssf_ar1(y), acf_linear=ssf_acf_linear(y),
                         spectral=ssf_spectral(y)))
    ep = pd.DataFrame(rows)
    part = ep.groupby("id")[["ar1", "acf_linear", "spectral"]].median()
    return dict(ar1=float(part.ar1.median()),
                acf_linear=float(part.acf_linear.median()),
                spectral=float(part.spectral.median()),
                n_participants_spectral=int(part.spectral.notna().sum()),
                n_intervals_spectral=int(ep.spectral.notna().sum()))


def instrument_ssf(df, data_dir):
    rows = []
    def add(label, d, col):
        rec = participant_balanced_ssf(d, col)
        rec["measure"] = label; rows.append(rec)

    add("E3G predictor", df.dropna(subset=["estrogen"]), "estrogen")
    for item in ("fatigue", "moodswing", "cramps", "bloating", "stress"):
        t = df[["id", "study_interval", "day_in_study", item]].copy()
        t["value"] = t[item].map(ORDINAL)
        add(item, t.dropna(subset=["value"]), "value")
    add("resting heart rate", daily_rhr(data_dir), "value")
    add("skin temperature", daily_temp(data_dir), "value")
    return pd.DataFrame(rows)


def differential_prediction(df):
    rows = []
    for item, cls in CLASSIFICATION.items():
        rec = item_records(df, item)
        pr = participant_r(rec)
        r_bal = fisher_aggregate(list(pr.values()))
        r_pool = pooled_episode_centered(rec)
        rows.append(dict(item=item, classification=cls,
                         participants=len(pr), intervals=len(rec),
                         participant_balanced_r=r_bal,
                         pooled_within_r=r_pool,
                         g_participant_balanced=effect_g(r_bal),
                         g_pooled=effect_g(r_pool)))
    return pd.DataFrame(rows)


def surrogate_test(df, item, rng, B=500):
    rec = item_records(df, item)
    obs = participant_r(rec)
    vals = np.asarray(list(obs.values()), float)
    S = float(np.std(vals))
    null = np.empty(B)
    for b in range(B):
        null[b] = np.std(list(participant_r(rec, rng).values()))
    return dict(item=item, participants=len(vals), intervals=len(rec),
                SD_obs=S, SD_null_median=float(np.median(null)),
                null_q025=float(np.quantile(null, .025)),
                null_q975=float(np.quantile(null, .975)),
                p=float((1 + np.sum(null >= S)) / (B + 1)))


def phase_locked(df):
    rows = []
    for item in ("fatigue", "moodswing", "cramps", "bloating", "sorebreasts"):
        d = df[["episode", "phase", "estrogen", item]].copy()
        d["y"] = d[item].map(ORDINAL)
        d = d.dropna(subset=["phase", "estrogen", "y"])
        d["ec"] = d.estrogen - d.groupby("episode").estrogen.transform("mean")
        def zfun(s):
            sd = s.std(ddof=0)
            return (s - s.mean()) / sd if sd > 1e-12 else s * 0.0
        d["yz"] = d.groupby("episode")["y"].transform(zfun)
        r_e3g = abs(corr0(d.ec, d.yz))
        grand = d.yz.mean()
        ss_total = float(((d.yz - grand) ** 2).sum())
        ss_between = float(sum(len(g) * (g.mean() - grand) ** 2 for _, g in d.groupby("phase").yz))
        eta2 = ss_between / ss_total if ss_total > 0 else np.nan
        mens = float(d[d.phase.astype(str).str.contains("enstrual", case=False, na=False)].yz.mean())
        rows.append(dict(item=item, abs_pooled_within_r_E3G=r_e3g,
                         phase_eta2_on_within_episode_z=eta2,
                         menstrual_mean_z=mens))
    return pd.DataFrame(rows)


def objective_records(df, data_dir, col):
    h = df[["id", "study_interval", "episode", "day_in_study", "estrogen", "fatigue"]].copy()
    h["fatigue_n"] = h.fatigue.map(ORDINAL)
    h = h.groupby(["id", "study_interval", "episode", "day_in_study"], as_index=False).agg(
        estrogen=("estrogen", "mean"), fatigue_n=("fatigue_n", "mean"))
    rhr = daily_rhr(data_dir).rename(columns={"value": "rhr"})
    temp = daily_temp(data_dir).rename(columns={"value": "temp"})
    m = (h.merge(rhr, on=["id", "study_interval", "day_in_study"], how="left")
          .merge(temp, on=["id", "study_interval", "day_in_study"], how="left"))
    out = []
    for (pid, interval, ep), g in m.dropna(subset=["estrogen", col]).groupby(["id", "study_interval", "episode"]):
        g = g.groupby("day_in_study", as_index=False).agg(x=("estrogen", "mean"), y=(col, "mean"))
        if len(g) < MIN_PAIRED:
            continue
        day = g.day_in_study.to_numpy(int); x = g.x.to_numpy(float); y = g.y.to_numpy(float)
        grid, pos = prepare_calendar(x, day)
        out.append(dict(id=int(pid), interval=interval, episode=ep,
                        x=x, y=y, day=day, n=len(g), grid=grid, pos=pos))
    return out


def objective_surrogate(df, data_dir, rng, B=500):
    rows = []
    for col, label in (("fatigue_n", "fatigue"), ("rhr", "resting heart rate"),
                       ("temp", "skin temperature")):
        rec = objective_records(df, data_dir, col)
        obs = participant_r(rec); vals = np.asarray(list(obs.values()), float)
        S = float(np.std(vals)); null = np.empty(B)
        for b in range(B):
            null[b] = np.std(list(participant_r(rec, rng).values()))
        rows.append(dict(outcome=label, participants=len(vals), intervals=len(rec),
                         SD_obs=S, SD_null_median=float(np.median(null)),
                         null_q025=float(np.quantile(null, .025)),
                         null_q975=float(np.quantile(null, .975)),
                         p=float((1 + np.sum(null >= S)) / (B + 1)),
                         median_abs_r=float(np.median(np.abs(vals)))))
    return pd.DataFrame(rows)


def attenuation_sensitivity(df, ssf, rng, n_boot=2000):
    sx = float(ssf.loc[ssf.measure == "E3G predictor", "spectral"].iloc[0])
    rows = []
    for item in ("fatigue", "moodswing", "cramps", "bloating"):
        rec = item_records(df, item); pr = participant_r(rec)
        vals = np.asarray(list(pr.values()), float)
        sy = float(ssf.loc[ssf.measure == item, "spectral"].iloc[0])
        r_obs = fisher_aggregate(vals)
        for f in (0.0, 0.25, 0.50, 0.75):
            rx = sx + f * (1 - sx); ry = sy + f * (1 - sy)
            att = np.sqrt(rx * ry); boot = np.empty(n_boot)
            for b in range(n_boot):
                rb = fisher_aggregate(vals[rng.integers(0, len(vals), len(vals))])
                boot[b] = effect_g(np.clip(rb / att, -0.99, 0.99))
            lo, hi = np.percentile(boot, [2.5, 97.5])
            rows.append(dict(item=item, f=f, ssf_predictor=sx, ssf_outcome=sy,
                             attenuation_sensitivity=att,
                             participant_balanced_r=r_obs,
                             g_sensitivity=effect_g(np.clip(r_obs / att, -0.99, 0.99)),
                             ci95_low=lo, ci95_high=hi))
    return pd.DataFrame(rows)


def legacy_objective_bug_audit(df, data_dir, rng, B=500):
    """Reproduce the legacy raw duplicate merge to quantify its effect. Audit only."""
    h = df.copy(); h["fatigue_n"] = h.fatigue.map(ORDINAL)
    rhr = pd.read_csv(os.path.join(data_dir, "resting_heart_rate.csv"))[
        ["id", "day_in_study", "value"]].rename(columns={"value": "rhr"})
    temp = pd.read_csv(os.path.join(data_dir, "computed_temperature.csv"))
    temp = temp[temp.type == "SKIN"][["id", "sleep_start_day_in_study", "nightly_temperature"]]
    temp = temp.rename(columns={"sleep_start_day_in_study": "day_in_study",
                                "nightly_temperature": "temp"})
    m = (h[["id", "day_in_study", "estrogen", "fatigue_n"]]
         .merge(rhr, on=["id", "day_in_study"], how="left")
         .merge(temp, on=["id", "day_in_study"], how="left"))
    rows = []
    for col in ("fatigue_n", "rhr", "temp"):
        per = []
        for _, g in m.groupby("id"):
            g = g.dropna(subset=["estrogen", col]).sort_values("day_in_study")
            if len(g) >= MIN_PAIRED:
                per.append((g.estrogen.to_numpy(float), g[col].to_numpy(float)))
        S = np.std([corr0(x, y) for x, y in per]); null = np.empty(B)
        for b in range(B):
            null[b] = np.std([corr0(phase_randomize(x, rng), y) for x, y in per])
        rows.append(dict(outcome=col, participants=len(per), SD_obs=S,
                         null_SD_median=float(np.median(null)),
                         p=float((1 + np.sum(null >= S)) / (B + 1)),
                         median_merged_rows_per_person=float(np.median([len(x) for x, y in per]))))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", default="./results")
    ap.add_argument("--seed", type=int, default=20260824)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--surrogates", type=int, default=500)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    df = load_hormones(a.data_dir)

    n_part = df.id.nunique()
    n_int = df[["id", "study_interval"]].drop_duplicates().shape[0]
    fat = item_records(df, "fatigue")
    print(f"mcPHASES: {n_part} participants, {n_int} participant-intervals; "
          f"{len(fat)} qualifying fatigue intervals; median paired={np.median([r['n'] for r in fat]):.0f}")

    ssf = instrument_ssf(df, a.data_dir)
    ssf.to_csv(os.path.join(a.out_dir, "instrument_ssf_participant_balanced.csv"), index=False)
    differential_prediction(df).to_csv(os.path.join(a.out_dir, "differential_episode_aware.csv"), index=False)
    pd.DataFrame([
        surrogate_test(df, "fatigue", rng, a.surrogates),
        surrogate_test(df, "moodswing", rng, a.surrogates),
    ]).to_csv(os.path.join(a.out_dir, "surrogate_confirmatory_episode_aware.csv"), index=False)
    phase_locked(df).to_csv(os.path.join(a.out_dir, "phase_locked_episode_z.csv"), index=False)
    objective_surrogate(df, a.data_dir, rng, a.surrogates).to_csv(
        os.path.join(a.out_dir, "objective_surrogate_episode_aware.csv"), index=False)
    attenuation_sensitivity(df, ssf, rng, a.boot).to_csv(
        os.path.join(a.out_dir, "attenuation_sensitivity_participant_cluster_bootstrap.csv"), index=False)
    legacy_objective_bug_audit(df, a.data_dir, rng, a.surrogates).to_csv(
        os.path.join(a.out_dir, "legacy_objective_bug_reproduction.csv"), index=False)

    print(ssf.to_string(index=False))
    print("[saved]", a.out_dir)


if __name__ == "__main__":
    main()
