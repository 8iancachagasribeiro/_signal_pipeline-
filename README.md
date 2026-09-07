# Reproducibility package: Article 1
## From Group Nulls to Person-Specific Recovery

This repository supports the manuscript **From Group Nulls to Person-Specific Recovery: Statistical Cancellation, Signal Availability, and Design in Intensive Longitudinal Research**.

## Canonical BRM release

The canonical submission-aligned record is the published GitHub Release:

https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.1

Use the release page, release asset, and release tag as the stable references for the BRM submission. Do **not** use a branch-tree URL as a permanent citation, because branches are mutable and this repository also contains a maintenance branch named `brm-v1.0.1`.

The validated archival simulation asset is:

`BRM_reproducibility_v1_v1.0.1.zip`

Its current SHA-256 checksum is recorded in `BRM_reproducibility_v1/RELEASE_ASSET_SHA256.txt` and in the release description.

## What the archival ZIP reproduces

The archival ZIP contains the complete BRM methodological-expansion package:

- a 400-cell robustness grid with 1,000 Monte Carlo replications per cell;
- a 384-cell sampling-design grid with 1,000 Monte Carlo replications per cell;
- an SSF benchmark based on 72 simulation cells and eight estimator variants, yielding 576 method-by-condition rows;
- a 31-cell phase-randomized surrogate calibration/power study with 1,000 Monte Carlo replications and 199 surrogates per test;
- deterministic seed maps, environment information, canonical result CSVs, summary tables, six manuscript figures, and validation utilities.

After extracting the release asset:

```bash
cd BRM_reproducibility_v1
python -m py_compile *.py
python validate_brm_outputs.py --root .
```

The validator checks canonical row counts, replication counts, and key manuscript numerical anchors.

## Empirical analysis code

The release tag also snapshots the audited empirical source code stored at repository level:

- `mcphases_analyses.py` — mcPHASES participant-interval processing, phase-randomized surrogate analyses, SSF summaries, predictor-alignment analyses, and sensitivity analyses;
- `actigraphy_replication.py` — independent clinical-actigraphy SSF replication with diagnosis labels read from source metadata;
- `ssf_estimators.py` — audited SSF estimators and explicit regular-grid gap handling.

The empirical scripts require the original source archives. Raw mcPHASES data are **not** redistributed because PhysioNet access is credentialed. Derived non-identifying audit outputs are retained under `audit_results/2026-08-24/` for provenance.

## OSF project and transparency

The public OSF project associated with the study is:

https://osf.io/4u6dk/

The OSF project provides project-level study materials and transparency information. It is not represented in the current manuscript as a formal preregistration.

## Effect-size notation

Correlations are converted for descriptive comparability using:

`d = 2r / sqrt(1 - r^2)`

The manuscript denotes these standardized values as **Cohen's d-equivalent effects**. No Hedges' g small-sample correction is applied.

## Recovery metrics

Recovery fidelity is a rank-ordering metric. The complementary ICC implemented in the canonical BRM simulations is a **single-measure consistency ICC, ICC(C,1)**, not an absolute-agreement ICC. RMSE and bias carry the magnitude-error interpretation.

## SSF is not reliability

The smooth-signal fraction (SSF) is a diagnostic of temporally structured signal availability at the target timescale. It is not classical reliability and does not by itself identify whether high-frequency variation is measurement error or genuine rapid biological variation.

## Data access

Raw mcPHASES data are not redistributed because PhysioNet access is credentialed. Independent actigraphy datasets remain available from their original repositories. See `BRM_reproducibility_v1/README_DATA.md` for the source records.

## Versioning

GitHub Release `brm-v1.0.1` is the corrected submission-aligned archival snapshot produced after the BRM link, metadata, and methodological-consistency audit. A permanent Zenodo DOI should be added only after an exact matching archival deposit is published and independently accessible.