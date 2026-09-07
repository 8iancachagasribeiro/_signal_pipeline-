#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def near(x, y, tol=5e-4):
    if not np.isfinite(x) or abs(float(x) - float(y)) > tol:
        raise AssertionError(f"{x} != {y} within {tol}")


def main():
    ap = argparse.ArgumentParser(description="Validate derived empirical audit anchors used by the BRM manuscript.")
    ap.add_argument("--root", default=".")
    a = ap.parse_args()
    root = Path(a.root).resolve()
    emp = root / "empirical_audit"

    focal = pd.read_csv(emp / "surrogate_focal_episode_aware.csv").set_index("item")
    obj = pd.read_csv(emp / "objective_surrogate_episode_aware.csv").set_index("outcome")
    ssf = pd.read_csv(emp / "instrument_ssf_participant_balanced.csv").set_index("measure")
    phase = pd.read_csv(emp / "phase_locked_episode_z.csv").set_index("item")

    assert int(focal.loc["fatigue", "participants"]) == 41
    assert int(focal.loc["moodswing", "participants"]) == 41
    near(focal.loc["fatigue", "SD_obs"], 0.1386594539, 1e-9)
    near(focal.loc["fatigue", "SD_null_median"], 0.1395427376, 1e-9)
    near(focal.loc["fatigue", "p"], 0.5329341317, 1e-9)
    near(focal.loc["moodswing", "SD_obs"], 0.1301228051, 1e-9)
    near(focal.loc["moodswing", "SD_null_median"], 0.1468251635, 1e-9)
    near(focal.loc["moodswing", "p"], 0.8802395210, 1e-9)

    near(obj.loc["resting heart rate", "SD_obs"], 0.2161018132, 1e-9)
    near(obj.loc["resting heart rate", "SD_null_median"], 0.1950408973, 1e-9)
    near(obj.loc["resting heart rate", "p"], 0.1556886228, 1e-9)
    near(obj.loc["skin temperature", "SD_obs"], 0.1155064279, 1e-9)
    near(obj.loc["skin temperature", "SD_null_median"], 0.1324407296, 1e-9)
    near(obj.loc["skin temperature", "p"], 0.8742514970, 1e-9)

    expected_ssf = {
        "E3G predictor": 0.4271236768,
        "fatigue": 0.3618242824,
        "moodswing": 0.3783188213,
        "cramps": 0.5724210897,
        "bloating": 0.3490580827,
        "resting heart rate": 0.8913937432,
        "skin temperature": 0.0819629845,
    }
    for name, value in expected_ssf.items():
        near(ssf.loc[name, "spectral"], value, 1e-9)

    near(phase.loc["cramps", "abs_pooled_within_r_E3G"], 0.0781559415, 1e-9)
    near(phase.loc["cramps", "phase_eta2_on_within_episode_z"], 0.1743926154, 1e-9)
    near(phase.loc["cramps", "menstrual_mean_z"], 0.8356445566, 1e-9)

    print("EMPIRICAL VALIDATION PASS: derived audit outputs match the manuscript empirical anchors.")


if __name__ == "__main__":
    main()
