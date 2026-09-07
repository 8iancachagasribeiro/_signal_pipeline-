# Release notes: v1.0.1-brm

This corrected version resolves transparency, repository-destination, versioning, and methodological-label inconsistencies identified during the BRM resubmission audit while preserving the validated numerical outputs.

## Included in the archival simulation asset

- audited reference model;
- 400-cell robustness grid, 1,000 Monte Carlo replications per cell;
- 384-cell sampling-design grid, 1,000 replications per cell;
- SSF benchmark across 72 simulation cells and eight estimator variants;
- surrogate-test calibration, homogeneous-effect boundary analysis, and power grid;
- deterministic seed map and exact execution environment;
- canonical result CSVs, summary tables, six generated figures, and validation utilities.

## Empirical source code in the release tag

The tagged repository source additionally includes the audited empirical scripts:

- `mcphases_analyses.py`;
- `actigraphy_replication.py`;
- `ssf_estimators.py`.

These scripts begin from the original source archives. Raw mcPHASES participant data are not redistributed because PhysioNet access is credentialed.

## Corrected destinations and metadata

- public OSF project: https://osf.io/4u6dk/;
- canonical release page: https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.1;
- branch-tree URLs are not used as permanent source citations because branches are mutable;
- correlation-derived standardized effects are described as Cohen's d-equivalent values rather than Hedges' g;
- no Zenodo DOI is claimed until a matching v1.0.1 archival deposit is published and independently accessible.

## Statistical-label correction

The complementary ICC produced by the canonical recovery code is a **single-measure consistency ICC, ICC(C,1)**. It is not an absolute-agreement ICC. Recovery fidelity remains an ordering metric, while RMSE and mean bias quantify magnitude error.

## Interpretation boundaries

- Recovery fidelity is an ordering metric, not proof of magnitude agreement.
- SSF is a smooth-signal availability diagnostic, not classical reliability.
- The phase-randomized surrogate procedure tests excess dispersion relative to its temporal surrogate null and is not a universal random-slope variance-component test under arbitrary nonzero homogeneous effects.

## Restricted data

mcPHASES raw data are not redistributed. See `README_DATA.md`.