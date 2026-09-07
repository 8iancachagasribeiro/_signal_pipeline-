# Changelog

## 2026-09-07 - BRM v1.0.2 final reproducibility hardening

### Canonical submission state
- Canonical public computational record: GitHub Release `brm-v1.0.2`.
- Public OSF project: https://osf.io/4u6dk/.
- Robustness grid: 400 cells, 1,000 Monte Carlo replications per cell.
- Sampling-design grid: 384 cells, 1,000 Monte Carlo replications per cell.
- SSF benchmark: 72 simulation cells x eight estimator variants = 576 method-by-condition rows, 1,000 replications per simulation cell.
- Surrogate calibration/power study: 31 cells, 1,000 Monte Carlo replications per cell, 199 surrogates per test.
- Current manuscript: six figures.

### Reproducibility corrections
- Added the audited mcPHASES and independent-actigraphy analysis scripts to the archival package.
- Added derived, non-identifying empirical audit outputs to permit direct numerical verification without redistributing restricted participant data.
- Added `validate_empirical_audit.py` for focal surrogate outcomes, objective outcomes, SSF anchors, and cramping phase-alignment anchors.
- Added `make_empirical_figure.py` to rebuild the empirical manuscript illustration from bundled derived audit outputs.
- Added explicit empirical-reproduction documentation, including the temporary interpolation used only to construct a regular grid for FFT phase randomization and the subsequent resampling only at originally observed days.
- Corrected the recovery-metric terminology: the implemented coefficient is `ICC(C,1)`, a single-measure consistency ICC, not an absolute-agreement ICC.
- Removed branch URLs as archival citations; the versioned release asset and its SHA-256 are the submission-facing computational record.
- Preserved the OSF destination as project-level transparency material without describing it as a formal preregistration.
- Retained correlation-derived standardized effects as Cohen's d-equivalent values with no Hedges' g correction.

### Numerical status
The validated numerical outputs are unchanged from the audited v1.0.1 computational results. Version 1.0.2 hardens packaging, terminology, empirical-code inclusion, and validation rather than changing the reported scientific results.

## 2026-09-05 - BRM v1.0.1 metadata/link correction
- Corrected the public OSF destination and removed inaccessible/incorrect Zenodo submission references.
- Expanded the principal Monte Carlo grids to the final 1,000-replication submission state.
- Published the validated v1.0.1 release asset and checksum.
- v1.0.1 is superseded by v1.0.2 for the current BRM resubmission.

## 2026-08-24 - full methodological and empirical audit
- Re-audited mcPHASES from the original credentialed archive with participant-interval temporal segmentation.
- Aggregated repeated intervals back to participant-level couplings by Fisher-z weighting rather than treating intervals as independent participants.
- Corrected daily Fitbit merging and removed a legacy duplicate-merge artifact affecting skin temperature.
- Confirmed null/focal surrogate results for fatigue and mood swing and null objective results for resting heart rate and skin temperature under the corrected pipeline.
- Re-estimated participant-balanced spectral SSF with strict gap handling.
- Retained cramping as the strongest phase-alignment illustration, with cycle phase carrying substantially more information than linear E3G level.
- Reframed SSF as a smooth-signal availability diagnostic rather than classical reliability.
- Added deterministic audit outputs and static/numerical consistency checks.
