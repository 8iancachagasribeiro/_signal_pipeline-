# Release notes: v1.0.2-brm

Version 1.0.2 is the final reproducibility-hardening release for the current BRM resubmission. It preserves the validated canonical numerical outputs while tightening terminology, empirical-code availability, versioning, and archival traceability.

## Included

- audited reference model;
- 400-cell robustness grid, 1,000 Monte Carlo replications per cell;
- 384-cell sampling-design grid, 1,000 replications per cell;
- SSF benchmark across 72 simulation cells and eight estimator variants;
- surrogate-test calibration, homogeneous-effect boundary analysis, and power grid;
- audited mcPHASES and actigraphy analysis scripts;
- derived non-identifying empirical audit outputs and empirical-anchor validator;
- deterministic seed map and exact execution environment;
- canonical result CSVs, summary tables, six generated figures, and validation utilities.

## Corrections and hardening in v1.0.2

- labels the reported ICC correctly as a single-measure consistency coefficient, `ICC(C,1)`, rather than an absolute-agreement ICC;
- includes the audited empirical mcPHASES and actigraphy pipelines in the archival package so Code Availability matches the actual release contents;
- documents temporary interpolation of internal missing predictor days solely for regular-grid FFT surrogate construction, followed by resampling only on originally observed days;
- removes reliance on a branch URL as the archival citation and uses the versioned GitHub Release as the canonical public computational record;
- bundles derived empirical audit outputs and validates the manuscript empirical numerical anchors;
- retains the public OSF project at https://osf.io/4u6dk/ as project-level transparency material without characterizing it as a formal preregistration;
- preserves Cohen's d-equivalent notation for correlation-derived standardized effects with no Hedges' g correction.

## Canonical destinations

- public OSF project: https://osf.io/4u6dk/;
- published v1.0.2 GitHub Release: https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.2;
- repository: https://github.com/8iancachagasribeiro/_signal_pipeline-;
- no Zenodo DOI is claimed until an exact matching v1.0.2 deposit is published and independently accessible.

## Interpretation boundaries

- Recovery fidelity is an ordering metric, not proof of magnitude agreement.
- `ICC(C,1)` measures consistency, not absolute agreement; RMSE and bias carry magnitude-error information.
- SSF is a smooth-signal availability diagnostic, not classical reliability.
- The phase-randomized surrogate procedure tests excess dispersion relative to its temporal surrogate null and is not a universal random-slope variance-component test under arbitrary nonzero homogeneous effects.

## Restricted data

mcPHASES raw data are not redistributed. See `README_DATA.md` and `EMPIRICAL_REPRODUCIBILITY.md`.