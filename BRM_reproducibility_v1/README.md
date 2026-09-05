# BRM reproducibility package v1.0.1

Repository companion for **From Group Nulls to Person-Specific Recovery: Statistical Cancellation, Signal Availability, and Design in Intensive Longitudinal Research**.

## Open practices

- Public analysis-plan materials: https://osf.io/4u6dk/
- Repository: https://github.com/8iancachagasribeiro/_signal_pipeline-
- Canonical v1.0.1 snapshot: https://github.com/8iancachagasribeiro/_signal_pipeline-/tree/brm-v1.0.1/BRM_reproducibility_v1
- Intended GitHub release tag: `brm-v1.0.1`
- Permanent Zenodo DOI: not assigned for v1.0.1

The manuscript distinguishes analyses specified in the initial analysis plan from subsequent simulation-based robustness and methodological validation analyses. The OSF project is not described here as a formal preregistration unless a separate timestamped, immutable OSF Registration is documented.

## Distribution model

This folder contains the executable BRM expansion code, metadata, compact manuscript tables, seeds, environment information, and validation utilities. The canonical submission destination is the versioned branch above. If an archival ZIP is published as a GitHub Release or deposited in Zenodo, it must be rebuilt from this v1.0.1 state and assigned a fresh checksum.

## Canonical analyses

1. `reference_model.py`: audited reference mechanism.
2. `brm_robustness.py`: 400 robustness cells, 1,000 Monte Carlo replications/cell.
3. `brm_sampling_design.py`: 384 sampling-design cells, 1,000 replications/cell.
4. `brm_ssf_benchmark.py`: 72 simulation cells and eight estimator variants, yielding 576 method rows.
5. `brm_surrogate_power.py`: zero-coupling calibration, homogeneous-effect boundary analysis, and 27 power cells; 1,000 replications and 199 surrogates/test.
6. `make_brm_outputs.py`: rebuilds summary tables plus PNG/SVG figures from canonical CSVs.
7. `validate_brm_outputs.py`: checks canonical row counts, replication counts, and manuscript numerical anchors.

## Validate the package

```bash
python -m py_compile *.py
python validate_brm_outputs.py --root .
```

After building any corrected archival ZIP, generate and verify a fresh SHA-256 before publishing it.

## Data access

Raw mcPHASES data are not redistributed because access is credentialed through PhysioNet. Independent actigraphy datasets remain at their original repositories. See `README_DATA.md`.

## Scientific boundaries

- Recovery fidelity measures preservation of individual ordering, not absolute magnitude agreement.
- SSF is a smooth-signal availability diagnostic, not classical reliability.
- The phase-randomized surrogate test detects excess coupling dispersion relative to the implemented temporal null and is not a universal random-slope variance-component test under arbitrary homogeneous nonzero effects.
- The Haar comparator is the specified implementation included here, not the full wavelet-method family.
- Correlation-derived standardized effects are reported as Cohen's d-equivalent values; no Hedges' g small-sample correction is applied.

## License

Repository-level license: GNU Affero General Public License v3.0.