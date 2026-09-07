# BRM reproducibility package v1.0.1

Repository companion for **From Group Nulls to Person-Specific Recovery: Statistical Cancellation, Signal Availability, and Design in Intensive Longitudinal Research**.

## Open practices

- Public OSF project: https://osf.io/4u6dk/
- Repository: https://github.com/8iancachagasribeiro/_signal_pipeline-
- Published GitHub Release: https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.1
- Permanent Zenodo DOI: not assigned for v1.0.1

The **GitHub Release page and its tag** are the canonical citation targets for this submission. A branch-tree URL is intentionally not used as a permanent source citation because branches are mutable.

The OSF project provides project-level study materials and transparency information. The executable reproducibility record for the current BRM submission is the versioned GitHub Release.

## Distribution model

This folder contains the executable BRM simulation expansion code, metadata, compact manuscript tables, seeds, environment information, and validation utilities. The published v1.0.1 GitHub Release contains the complete validated archival ZIP, including full canonical outputs and figures. The release checksum is recorded in `RELEASE_ASSET_SHA256.txt`.

The **tagged repository source** additionally contains the audited empirical scripts `mcphases_analyses.py`, `actigraphy_replication.py`, and `ssf_estimators.py`. Those scripts require the original source archives; restricted mcPHASES participant data are not included in the archival ZIP or repository.

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

For the published archival ZIP, also validate the internal checksum manifest after extraction.

## Data access

Raw mcPHASES data are not redistributed because access is credentialed through PhysioNet. Independent actigraphy datasets remain at their original repositories. See `README_DATA.md`.

## Scientific boundaries

- Recovery fidelity measures preservation of individual ordering, not absolute magnitude agreement.
- The ICC implemented by `rep_metrics()` is a **single-measure consistency ICC, ICC(C,1)**. It should not be described as an absolute-agreement ICC.
- RMSE and mean bias quantify magnitude error and systematic scale/location error.
- SSF is a smooth-signal availability diagnostic, not classical reliability.
- The phase-randomized surrogate test detects excess coupling dispersion relative to the implemented temporal null and is not a universal random-slope variance-component test under arbitrary homogeneous nonzero effects.
- The Haar comparator is the specified implementation included here, not the full wavelet-method family.
- Correlation-derived standardized effects are reported as Cohen's d-equivalent values; no Hedges' g small-sample correction is applied.

## License

Repository-level license: GNU Affero General Public License v3.0.