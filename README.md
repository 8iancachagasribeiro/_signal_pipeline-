# Reproducibility package: Article 1
## From Group Nulls to Person-Specific Recovery

This repository supports the manuscript **From Group Nulls to Person-Specific Recovery: Statistical Cancellation, Signal Availability, and Design in Intensive Longitudinal Research**.

## Canonical BRM reproduction package

The submission-aligned reproducibility package is located in `BRM_reproducibility_v1/`. The canonical versioned snapshot for the current BRM submission is:

https://github.com/8iancachagasribeiro/_signal_pipeline-/tree/brm-v1.0.1/BRM_reproducibility_v1

That folder, rather than older root-level development scripts or files under `archive/`, is the canonical source for the BRM methodological expansion.

The canonical package contains:

- a 400-cell robustness grid with 1,000 Monte Carlo replications per cell;
- a 384-cell sampling-design grid with 1,000 replications per cell;
- an SSF benchmark based on 72 simulation cells and eight estimator variants, yielding 576 method-by-condition rows;
- a 31-cell phase-randomized surrogate calibration/power study with 1,000 Monte Carlo replications and 199 surrogates per test;
- deterministic seed maps, environment information, summary tables, six manuscript figures, and validation utilities.

Run the package validator with:

```bash
cd BRM_reproducibility_v1
python validate_brm_outputs.py --root .
```

The validator checks the canonical row counts, replication counts, and key manuscript numerical anchors.

## Analysis-plan materials and transparency

The public OSF project containing the analysis-plan materials is:

https://osf.io/4u6dk/

The manuscript distinguishes analyses specified in the initial analysis plan from subsequent simulation-based robustness and methodological validation analyses. This repository does not characterize the OSF project itself as a formal immutable preregistration unless a separate timestamped OSF Registration is documented.

## Effect-size notation

Correlations are converted for descriptive comparability using:

`d = 2r / sqrt(1 - r^2)`

The manuscript denotes these standardized values as **Cohen's d-equivalent effects**. No Hedges' g small-sample correction is applied.

## SSF is not reliability

The smooth-signal fraction (SSF) is a diagnostic of temporally structured signal availability at the target timescale. It is not classical reliability and does not by itself identify whether high-frequency variation is measurement error or genuine rapid biological variation.

## Data access

Raw mcPHASES data are not redistributed because PhysioNet access is credentialed. Independent actigraphy datasets remain available from their original repositories. See `BRM_reproducibility_v1/README_DATA.md` for the source records.

## Versioning

`brm-v1.0.1` is the corrected submission-aligned branch produced after the BRM link and metadata audit. The historical `brm-v1.0.0` release is superseded for submission purposes because its metadata contained obsolete repository destinations. A permanent Zenodo DOI should be added only after a matching v1.0.1 archival deposit is published and independently accessible.

If a GitHub Release is created for this version, its intended tag is `brm-v1.0.1`.