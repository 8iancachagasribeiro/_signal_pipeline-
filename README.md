# Reproducibility package: Article 1
## From Group Nulls to Person-Specific Recovery

This repository supports the manuscript **From Group Nulls to Person-Specific Recovery: Statistical Cancellation, Signal Availability, and Design in Intensive Longitudinal Research**.

## Canonical BRM release

The submission-aligned archival release is:

https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.2

Release asset: `BRM_reproducibility_v1_v1.0.2.zip`

The release asset, rather than mutable development files elsewhere in the repository, is the canonical computational record for the current BRM submission. Its SHA-256 is recorded in `BRM_reproducibility_v1/RELEASE_ASSET_SHA256.txt` and in the GitHub Release metadata.

The canonical package contains:

- a 400-cell robustness grid with 1,000 Monte Carlo replications per cell;
- a 384-cell sampling-design grid with 1,000 Monte Carlo replications per cell;
- an SSF benchmark based on 72 simulation cells and eight estimator variants, yielding 576 method-by-condition rows;
- a 31-cell phase-randomized surrogate calibration/power study with 1,000 Monte Carlo replications and 199 surrogates per test;
- audited mcPHASES and actigraphy analysis scripts for reproduction from the original source archives;
- derived non-identifying empirical audit outputs used to verify manuscript anchors;
- deterministic seed maps, environment information, summary tables, six manuscript figures, and validation utilities.

After extracting the release asset:

```bash
cd BRM_reproducibility_v1
python -m py_compile *.py
python validate_brm_outputs.py --root .
python validate_empirical_audit.py --root .
sha256sum -c SHA256SUMS.txt
```

## OSF project and transparency

The public OSF project associated with the study is:

https://osf.io/4u6dk/

The OSF project provides project-level study materials and transparency information. It is not described here as a formal preregistration.

## Effect-size notation

Correlations are converted for descriptive comparability using:

`d = 2r / sqrt(1 - r^2)`

The manuscript denotes these standardized values as **Cohen's d-equivalent effects**. No Hedges' g small-sample correction is applied.

## Recovery metrics

Recovery fidelity is a correlation-based ordering metric. The package also reports RMSE, mean bias, directional accuracy, and a single-measure consistency ICC, `ICC(C,1)`. The ICC is not an absolute-agreement coefficient; magnitude recovery is evaluated jointly with RMSE and bias.

## SSF is not reliability

The smooth-signal fraction (SSF) is a diagnostic of temporally structured signal availability at the target timescale. It is not classical reliability and does not by itself identify whether high-frequency variation is measurement error or genuine rapid biological variation.

## Data access

Raw mcPHASES data are not redistributed because PhysioNet access is credentialed. Independent actigraphy datasets remain available from their original repositories. See `BRM_reproducibility_v1/README_DATA.md` and `BRM_reproducibility_v1/EMPIRICAL_REPRODUCIBILITY.md`.

## Versioning

GitHub Release `brm-v1.0.2` supersedes `brm-v1.0.1` for the BRM resubmission. Version 1.0.2 hardens empirical-code inclusion, terminology, validation, and archival citation while preserving the validated canonical numerical outputs. A Zenodo DOI should be cited only after an exact matching deposit is publicly available.