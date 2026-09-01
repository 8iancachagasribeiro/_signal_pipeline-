# Release notes: v1.0.0-brm

This snapshot packages the methodological expansion prepared for submission to *Behavior Research Methods*.

## Included

- audited reference model;
- 400-cell robustness grid, 1,000 Monte Carlo replications per cell;
- 384-cell sampling-design grid, 1,000 replications per cell;
- SSF benchmark across 72 simulation cells and eight estimator variants;
- surrogate-test calibration, homogeneous-effect boundary analysis, and power grid;
- deterministic seed map and exact execution environment;
- canonical result CSVs, summary tables, six generated figures, and validation utilities.

## Interpretation boundaries

- Recovery fidelity is an ordering metric, not proof of magnitude agreement.
- SSF is a smooth-signal availability diagnostic, not classical reliability.
- The phase-randomized surrogate procedure tests excess dispersion relative to its temporal surrogate null and is not a universal random-slope variance-component test under arbitrary nonzero homogeneous effects.
- Expanded robustness analyses are post-preregistration methodological validation analyses. The preregistration is at https://osf.io/ewyp7.

## Restricted data

mcPHASES raw data are not redistributed. See `README_DATA.md`.
