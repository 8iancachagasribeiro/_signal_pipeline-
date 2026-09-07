# Empirical reproducibility

This document describes the empirical pipelines bundled with BRM reproducibility package v1.0.2. The scripts are included so that credentialed users can begin from the original source archives rather than from manuscript numbers.

## mcPHASES

Raw mcPHASES data are not redistributed. Obtain the dataset from PhysioNet under its credentialed-access conditions (DOI: https://doi.org/10.13026/zx6a-2c81).

The audited script expects a source directory containing the files used by the published archive, including:

- `hormones_and_selfreport.csv`
- `resting_heart_rate.csv`
- `computed_temperature.csv`

Run:

```bash
python mcphases_analyses.py --data-dir /path/to/mcphases --out-dir empirical_results --seed 20260824 --surrogates 500 --boot 2000
```

Audit rules implemented in the script include:

- each `study_interval` is treated as a separate temporal segment;
- intervals separated by the long calendar gap are never bridged for FFT/surrogate operations;
- repeated intervals are precision-weighted with Fisher transformation back to one participant-level coupling estimate;
- daily Fitbit measures are collapsed to one participant-interval-day before merging;
- resting-heart-rate values <= 0 are treated as missing;
- internal missing predictor days are linearly interpolated only to construct the regular daily grid required for FFT phase randomization;
- each generated surrogate is sampled back only on originally observed days;
- SSF is treated as a smooth-signal availability diagnostic, not classical reliability.

The release also bundles derived, non-identifying audit outputs under `empirical_audit/`. These permit numerical verification of manuscript anchors without redistributing raw participant data.

## Independent actigraphy

The source datasets remain at their original repositories:

- DEPRESJON: https://doi.org/10.5281/zenodo.1219550
- PSYKOSE: https://osf.io/dgjzu
- HYPERAKTIV: https://osf.io/3agwr

After arranging the source directories as described by their repositories, run:

```bash
python actigraphy_replication.py --data-dir /path/to/actigraphy_sources --out-dir actigraphy_results
```

HYPERAKTIV diagnostic labels are read from `patient_info.csv`; unresolved IDs are excluded rather than inferred from folder membership.

## Empirical figure and anchor validation

The bundled derived outputs can be checked without source data:

```bash
python validate_empirical_audit.py --root .
python make_empirical_figure.py --root . --out figures/figure6_empirical_illustration_rebuilt.png
```

The validator checks the focal fatigue and mood-swing surrogate results, objective outcomes, participant-balanced spectral SSF values, and the cramping phase-alignment anchors used in the manuscript.

## Scope

Successful execution depends on the source archives retaining the variable/file structure documented by their published versions. No script can bypass PhysioNet credential requirements, and no restricted data are included in this package.