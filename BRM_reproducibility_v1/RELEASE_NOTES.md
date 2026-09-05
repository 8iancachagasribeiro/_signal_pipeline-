# Release notes: v1.0.1-brm

This corrected version resolves transparency, repository-destination, and versioning inconsistencies identified during the BRM resubmission audit while preserving the validated numerical outputs.

## Included

- audited reference model;
- 400-cell robustness grid, 1,000 Monte Carlo replications per cell;
- 384-cell sampling-design grid, 1,000 replications per cell;
- SSF benchmark across 72 simulation cells and eight estimator variants;
- surrogate-test calibration, homogeneous-effect boundary analysis, and power grid;
- deterministic seed map and exact execution environment;
- canonical result CSVs, summary tables, six generated figures, and validation utilities.

## Corrected destinations and metadata

- public analysis-plan materials: https://osf.io/4u6dk/;
- canonical v1.0.1 snapshot: https://github.com/8iancachagasribeiro/_signal_pipeline-/tree/brm-v1.0.1/BRM_reproducibility_v1;
- the OSF project is not labeled as a formal preregistration unless a separate immutable Registration is documented;
- correlation-derived standardized effects are described as Cohen's d-equivalent values rather than Hedges' g;
- the package documents 1,000 Monte Carlo replications for the canonical BRM grids and six manuscript figures;
- no Zenodo DOI is claimed until a matching v1.0.1 archival deposit is published and independently accessible.

## Interpretation boundaries

- Recovery fidelity is an ordering metric, not proof of magnitude agreement.
- SSF is a smooth-signal availability diagnostic, not classical reliability.
- The phase-randomized surrogate procedure tests excess dispersion relative to its temporal surrogate null and is not a universal random-slope variance-component test under arbitrary nonzero homogeneous effects.

## Restricted data

mcPHASES raw data are not redistributed. See `README_DATA.md`.