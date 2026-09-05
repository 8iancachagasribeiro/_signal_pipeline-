# Reproducibility package: Article 1
## From Group Nulls to Person-Specific Recovery

This repository supports the manuscript **From Group Nulls to Person-Specific Recovery: Statistical Cancellation, Signal Availability, and Design in Intensive Longitudinal Research**.

## Canonical BRM reproduction package

The corrected submission-aligned GitHub Release for the current BRM manuscript is:

https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.1

The corresponding source snapshot is located in `BRM_reproducibility_v1/` on branch `brm-v1.0.1`:

https://github.com/8iancachagasribeiro/_signal_pipeline-/tree/brm-v1.0.1/BRM_reproducibility_v1

These versioned materials, rather than older root-level development scripts or files under `archive/`, are the canonical source for the BRM methodological expansion.

The canonical package contains:

- a 400-cell robustness grid with 1,000 Monte Carlo replications per cell;
- a 384-cell sampling-design grid with 1,000 Monte Carlo replications per cell;
- an SSF benchmark based on 72 simulation cells and eight estimator variants, yielding 576 method-by-condition rows;
- a 31-cell phase-randomized surrogate calibration/power study with 1,000 Monte Carlo replications and 199 surrogates per test;
- deterministic seed maps, environment information, summary tables, six manuscript figures, and validation utilities.

Release asset: `BRM_reproducibility_v1_v1.0.1.zip`

The current release checksum is recorded in `BRM_reproducibility_v1/RELEASE_ASSET_SHA256.txt`.

Run the package validator after extraction with:

```bash
cd BRM_reproducibility_v1
python validate_brm_outputs.py --root .
```

The validator checks the canonical row counts, replication counts, and key manuscript numerical anchors.

## OSF project and transparency

The public OSF project associated with the study is:

https://osf.io/4u6dk/

The OSF project provides project-level study materials and transparency information. The canonical executable reproducibility materials are the versioned GitHub Release and source snapshot linked above.

## Effect-size notation

Correlations are converted for descriptive comparability using:

`d = 2r / sqrt(1 - r^2)`

The manuscript denotes these standardized values as **Cohen's d-equivalent effects**. No Hedges' g small-sample correction is applied.

## SSF is not reliability

The smooth-signal fraction (SSF) is a diagnostic of temporally structured signal availability at the target timescale. It is not classical reliability and does not by itself identify whether high-frequency variation is measurement error or genuine rapid biological variation.

## Data access

Raw mcPHASES data are not redistributed because PhysioNet access is credentialed. Independent actigraphy datasets remain available from their original repositories. See `BRM_reproducibility_v1/README_DATA.md` for the source records.

## Versioning

GitHub Release `brm-v1.0.1` is the corrected submission-aligned archival snapshot produced after the BRM link and metadata audit. The historical `brm-v1.0.0` release is superseded for submission purposes because its metadata contained obsolete repository destinations. A permanent Zenodo DOI should be added only after a matching v1.0.1 archival deposit is published and independently accessible.