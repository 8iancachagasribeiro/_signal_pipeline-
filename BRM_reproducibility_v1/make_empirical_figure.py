#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(description="Rebuild the BRM empirical illustration from derived non-identifying audit outputs.")
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default="figures/figure6_empirical_illustration_rebuilt.png")
    a = ap.parse_args()
    root = Path(a.root).resolve()
    emp_dir = root / "empirical_audit"

    focal = pd.read_csv(emp_dir / "surrogate_focal_episode_aware.csv")
    obj = pd.read_csv(emp_dir / "objective_surrogate_episode_aware.csv")
    ssf = pd.read_csv(emp_dir / "instrument_ssf_participant_balanced.csv")
    phase = pd.read_csv(emp_dir / "phase_locked_episode_z.csv")

    rows = []
    for label, key in (("fatigue", "fatigue"), ("mood swing", "moodswing")):
        r = focal.loc[focal.item == key].iloc[0]
        rows.append((label, r.SD_obs, r.SD_null_median, r.p))
    for label in ("resting heart rate", "skin temperature"):
        r = obj.loc[obj.outcome == label].iloc[0]
        rows.append((label, r.SD_obs, r.SD_null_median, r.p))
    het = pd.DataFrame(rows, columns=["outcome", "observed_sd", "null_median", "p"])

    measure_order = ["E3G predictor", "fatigue", "moodswing", "cramps", "bloating", "resting heart rate", "skin temperature"]
    ssf_plot = ssf.set_index("measure").loc[measure_order].reset_index()
    ssf_plot["display"] = ssf_plot.measure.replace({"E3G predictor": "E3G", "moodswing": "mood swing"})

    cr = phase.loc[phase.item == "cramps"].iloc[0]

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    xx = np.arange(len(het)); w = .36
    ax[0].bar(xx - w/2, het.observed_sd, w, label="Observed SD(r_i)")
    ax[0].bar(xx + w/2, het.null_median, w, label="Null median")
    ax[0].set_xticks(xx); ax[0].set_xticklabels(het.outcome, rotation=25, ha="right")
    ax[0].set_ylabel("Between-person dispersion")
    ax[0].set_title("(a) Coupling heterogeneity")
    ax[0].legend(fontsize=8)
    for i, p in enumerate(het.p):
        ax[0].text(i, max(het.observed_sd.iloc[i], het.null_median.iloc[i]) + .006, f"p={p:.3f}", ha="center", fontsize=7)

    xx = np.arange(len(ssf_plot))
    ax[1].bar(xx, ssf_plot.spectral)
    ax[1].set_xticks(xx); ax[1].set_xticklabels(ssf_plot.display, rotation=30, ha="right", fontsize=8)
    ax[1].set_ylim(0, 1); ax[1].set_ylabel("SSF"); ax[1].set_title("(b) Empirical smooth-signal fraction")

    ax[2].bar([0, 1], [cr.abs_pooled_within_r_E3G, cr.phase_eta2_on_within_episode_z])
    ax[2].set_xticks([0, 1]); ax[2].set_xticklabels(["|r| with E3G level", "Phase eta-squared"])
    ax[2].set_ylim(0, .22); ax[2].set_title("(c) Cramps: predictor alignment")

    fig.tight_layout()
    out = root / a.out
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220, bbox_inches="tight")
    svg = out.with_suffix(".svg")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out}")
    print(f"[saved] {svg}")


if __name__ == "__main__":
    main()
