# Reproducibility package: Article 1
## Group nulls, intra-individual heterogeneity, and measurement limits

This repository contains the simulation and empirical-analysis pipeline supporting Article 1. The audit separates three quantities that must not be conflated: **classical reliability**, **smooth-signal fraction (SSF)**, and **the fraction of smooth outcome variance actually coupled to E2**.

The third quantity is not identified by SSF alone. The pipeline therefore no longer reports one value called “power with the real instruments.” Power based on empirical SSF is a **sensitivity analysis** over `q`, the share of smooth outcome variance attributable to the E2-coupled mechanism.

## Repository scope

The root directory contains the **canonical manuscript reproduction pipeline**. Derived outputs are kept under `results/` and dated audit outputs under `audit_results/`. Exploratory or superseded development scripts are retained under `archive/` for provenance and are **not part of the canonical reproduction order** unless explicitly referenced below.

## Environment

Tested with Python 3.12.

```bash
pip install -r requirements.txt
```

## Reproduction order

```bash
# 0. Package integrity and smoke tests
python audit_consistency.py

# 1. Validate SSF estimators
python ssf_estimators.py

# 2. Reproduce masking mechanism and directional-imbalance boundary
python h4_frontier.py --validate-only

# 3. Null calibration, recovery fidelity, aliasing, 500 replicates/cell
python calibration_fidelity_aliasing.py --out-dir ./results

# 4. Full theoretical H4 grid at 500 replicates
bash rodar_grades_500.sh

# 5. SSF-calibrated registered-test sensitivity
python registered_test_power.py diagnostic
python registered_test_power.py calib --out-dir ./results
python registered_test_power.py power --out-dir ./results

# 6. Attenuation-identification and budget sensitivity
python budget_allocation.py --table 9 --out-dir ./results
python budget_allocation.py --table 17 --out-dir ./results
python h4_v2.py grid --out-dir ./results
python h4_v2.py budget --out-dir ./results

# 7. Alternative q x heterogeneity sensitivity
python sigma_sweep.py --out-dir ./results

# 8. EMPIRICAL mcPHASES analyses, credentialed PhysioNet data required
python mcphases_analyses.py --data-dir /path/to/mcphases --out-dir ./results

# 9. Independent clinical-actigraphy replication
python actigraphy_replication.py --data-dir /path/to/actigraphy --out-dir ./results

# 10. Regenerate the five current manuscript figures
python make_figures.py --results ./results --out ./figures
```

## Script map

| Script | Role |
|---|---|
| `h4_frontier.py` | Generative estradiol → dopaminergic tone → inverted-U performance mechanism; masking and directional-imbalance boundary. |
| `calibration_fidelity_aliasing.py` | Theoretical false-positive behavior, recovery fidelity, and phase-targeted vs evenly spaced sampling. Uses **reliability** as a theoretical design parameter, not empirical SSF. |
| `ssf_estimators.py` | AR(1), ACF-linear and adopted spectral SSF estimators; validation against known smooth-signal fractions; gap-aware utilities. |
| `ssf_power.py` | Direct SSF calibration utilities; prevents the former SSF-as-reliability error. |
| `registered_test_power.py` | Phase-randomized surrogate-test sensitivity under empirical SSF targets and explicit `q`; regular calendar before masking observation days. |
| `fastlrt.py` | Random-slope likelihood-ratio tests for common and participant-specific predictor trajectories. |
| `h4_v2.py` | Predictor/outcome SSF and research-budget sensitivity. No single empirical-power interpretation. |
| `budget_allocation.py` | Classical attenuation assumption analysis over `f` and design sensitivity over `q`. |
| `mcphases_analyses.py` | Empirical mcPHASES analyses; person-centered group summaries, SSF attenuation sensitivity, calendar-aware phase randomization. |
| `actigraphy_replication.py` | Clinical-actigraphy SSF replication with corrected HYPERAKTIV labels, de-duplication, gap-free segments, and saved outputs. |
| `make_figures.py` | Implements all five current Article 1 figures and skips rather than fabricates missing empirical figures. |
| `audit_consistency.py` | Static and numerical smoke-test audit. |

## Core reproducible claim

The generative mechanism remains reproducible after audit. With the calibrated mechanism around `sigma_b = 0.05`, the validation run returns a group effect near `|g| = 0.081` while the typical individual effect is near `|g| = 0.32`, with positive and negative individual slopes in both directions.

```bash
python h4_frontier.py --validate-only
```

This mechanism does **not** depend on the later SSF-based power analyses.

## SSF is not reliability

The spectral SSF estimates how much observed variance is carried by temporal structure below the specified frequency cutoff. It classifies high-frequency variance as non-smooth regardless of whether that variance reflects measurement error or genuine rapid biological or affective variation.

```text
SSF != classical reliability
```

The classical disattenuation expression `sqrt(Rx * Ry)` can be explored only after an explicit reliability assumption. `budget_allocation.py --table 9` parameterizes this by `f`, the fraction of high-frequency variance treated as genuine process.

## Why power is a sensitivity surface

Even if outcome SSF is known exactly, it does not reveal what fraction of its smooth variance is driven by E2. A smooth outcome may contain other low-frequency processes. The simulator therefore defines:

```text
q = share of smooth outcome variance attributable to the E2-coupled mechanism
```

`registered_test_power.py` and `h4_v2.py` vary `q` instead of pretending that SSF alone identifies empirical power.

## mcPHASES gap handling

Phase randomization requires a regular time axis. Earlier code compacted missing paired days before the FFT. The audited code places the predictor on its integer-day calendar first, fills gaps only to construct the surrogate spectrum, randomizes phases on that regular grid, and samples the surrogate back at the originally observed days. Outcome observations are never invented.

The corrected empirical analyses were regenerated from the credentialed mcPHASES archive during the 2026-08-24 audit. The derived non-identifying results and methodological consequences are documented in `audit_results/2026-08-24/` and `CHANGELOG.md`.

## Actigraphy replication

HYPERAKTIV's activity folder mixes ADHD participants and controls. The current script reads labels from `patient_info.csv`; unresolved subjects are excluded rather than guessed. The script saves the exact analysis table and cycle-count robustness table, so the final N must be read from the current output rather than copied from an older manuscript version.

## Figure numbering

The current builder matches the audited Article 1 structure:

1. estimator validation, empirical SSF, objective-vs-self-report check;
2. masking mechanism, recovery fidelity, sparse-sample false-positive behavior;
3. attenuation-identification sensitivity and SSF-calibrated power sensitivity;
4. hormone level versus cycle phase;
5. instrument/design sensitivity over `q`.

## Data availability

mcPHASES requires credentialed PhysioNet access and is not redistributed here. The actigraphy datasets are externally available from their original repositories. See `README_DATA.txt` and the manuscript for source citations.

## Repository and archive versioning

This GitHub repository is the canonical development version. Add a release DOI to the manuscript only after the audited branch is merged and the matching snapshot is deposited. Do not cite an older Zenodo DOI for a different code state.
