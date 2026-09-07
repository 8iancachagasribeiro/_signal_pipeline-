# BRM reproducibility package v1.0.2

Repository companion for **From Group Nulls to Person-Specific Recovery: Statistical Cancellation, Signal Availability, and Design in Intensive Longitudinal Research**.

## Open practices

- Public OSF project: https://osf.io/4u6dk/
- Repository: https://github.com/8iancachagasribeiro/_signal_pipeline-
- Published GitHub Release: https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.2
- Permanent Zenodo DOI: not assigned for v1.0.2

The OSF project provides project-level study materials and transparency information. The executable reproducibility record for the current BRM submission is the versioned GitHub Release asset above.

## Distribution model

The published v1.0.2 release contains the complete validated archival ZIP, including canonical Monte Carlo outputs, six figures, empirical audit outputs, empirical analysis scripts, metadata, seeds, environment information, compact manuscript tables, and validation utilities. The outer archive checksum is recorded in `RELEASE_ASSET_SHA256.txt` outside the archive and in the GitHub Release metadata; `SHA256SUMS.txt` inside the archive validates its packaged files.

## Canonical simulation analyses

1. `reference_model.py`: audited reference mechanism.
2. `brm_robustness.py`: 400 robustness cells, 1,000 Monte Carlo replications/cell.
3. `brm_sampling_design.py`: 384 sampling-design cells, 1,000 replications/cell.
4. `brm_ssf_benchmark.py`: 72 simulation cells and eight estimator variants, yielding 576 method rows.
5. `brm_surrogate_power.py`: zero-coupling calibration, homogeneous-effect boundary analysis, and 27 power cells; 1,000 replications and 199 surrogates/test.
6. `make_brm_outputs.py`: rebuilds summary tables plus simulation/summary figures from canonical CSVs.
7. `validate_brm_outputs.py`: checks canonical row counts, replication counts, and manuscript numerical anchors.

## Empirical reproduction

- `mcphases_analyses.py`: audited mcPHASES analysis pipeline beginning from the original credentialed PhysioNet archive.
- `actigraphy_replication.py`: independent actigraphy SSF replication beginning from the source repositories.
- `ssf_estimators.py`: SSF estimators and explicit gap-handling utilities used by the empirical pipelines.
- `make_empirical_figure.py`: rebuilds the manuscript empirical illustration from the derived, non-identifying audit outputs bundled under `empirical_audit/`.
- `validate_empirical_audit.py`: verifies the empirical numerical anchors reported in the manuscript against the bundled derived audit outputs.

See `EMPIRICAL_REPRODUCIBILITY.md` for required source files and commands. Raw mcPHASES data are not redistributed.

## Validate the package

```bash
python -m py_compile *.py
python validate_brm_outputs.py --root .
python validate_empirical_audit.py --root .
sha256sum -c SHA256SUMS.txt
```

## Data access

Raw mcPHASES data are not redistributed because access is credentialed through PhysioNet. Independent actigraphy datasets remain at their original repositories. See `README_DATA.md`.

## Scientific boundaries

- Recovery fidelity measures preservation of individual ordering, not absolute magnitude agreement.
- The reported ICC is a single-measure consistency coefficient, `ICC(C,1)`, not an absolute-agreement ICC.
- SSF is a smooth-signal availability diagnostic, not classical reliability.
- The phase-randomized surrogate test detects excess coupling dispersion relative to the implemented temporal null and is not a universal random-slope variance-component test under arbitrary homogeneous nonzero effects.
- The Haar comparator is the specified implementation included here, not the full wavelet-method family.
- Correlation-derived standardized effects are reported as Cohen's d-equivalent values; no Hedges' g small-sample correction is applied.

## License

Repository-level license: GNU Affero General Public License v3.0.