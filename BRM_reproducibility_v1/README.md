# BRM reproducibility package v1.0.0

Release-ready repository companion for **From Group Nulls to Person-Specific Recovery: Statistical Cancellation, Signal Availability, and Design in Intensive Longitudinal Research**.

## Open practices

- Preregistration: https://osf.io/ewyp7
- Repository: https://github.com/8iancachagasribeiro/_signal_pipeline-
- BRM branch: https://github.com/8iancachagasribeiro/_signal_pipeline-/tree/brm-methodological-expansion
- Planned release tag: `brm-v1.0.0`
- Permanent Zenodo DOI: to be assigned after archiving the GitHub release.

The expanded robustness and validation analyses are post-preregistration methodological analyses and are explicitly distinguished from preregistered analyses.

## Distribution model

This GitHub folder contains the executable BRM expansion code, metadata, compact manuscript tables, seeds, environment information, and validation utilities. The full canonical Monte Carlo CSVs and publication-ready figures are distributed in the validated release asset **`BRM_reproducibility_v1.zip`**, whose SHA-256 is recorded in `RELEASE_ASSET_SHA256.txt`.

The release asset is the immutable archival snapshot to attach to the GitHub Release and deposit in Zenodo. Keeping the large derived outputs in the release asset avoids duplicating generated data in Git history while preserving a hash-verifiable research record.

## Canonical analyses

1. `reference_model.py`: audited reference mechanism.
2. `brm_robustness.py`: 400 robustness cells, 1,000 Monte Carlo replications/cell.
3. `brm_sampling_design.py`: 384 sampling-design cells, 1,000 replications/cell.
4. `brm_ssf_benchmark.py`: 72 simulation cells and eight estimator variants, yielding 576 method rows.
5. `brm_surrogate_power.py`: zero-coupling calibration, homogeneous-effect boundary analysis, and 27 power cells; 1,000 replications and 199 surrogates/test.
6. `make_brm_outputs.py`: rebuilds summary tables plus PNG/SVG figures from canonical CSVs.
7. `validate_brm_outputs.py`: checks canonical row counts, replication counts, and manuscript numerical anchors.

## Validate the release asset

```bash
sha256sum -c RELEASE_ASSET_SHA256.txt
unzip BRM_reproducibility_v1.zip
cd BRM_reproducibility_v1
python -m py_compile *.py
python validate_brm_outputs.py --root .
sha256sum -c SHA256SUMS.txt
```

The final local release archive passed all four checks before synchronization.

## Data access

Raw mcPHASES data are not redistributed because access is credentialed through PhysioNet. Independent actigraphy datasets remain at their original repositories. See `README_DATA.md`.

## Scientific boundaries

- Recovery fidelity measures preservation of individual ordering, not absolute magnitude agreement.
- SSF is a smooth-signal availability diagnostic, not classical reliability.
- The phase-randomized surrogate test detects excess coupling dispersion relative to the implemented temporal null and is not a universal random-slope variance-component test under arbitrary homogeneous nonzero effects.
- The Haar comparator is the specified NumPy implementation included here, not the full wavelet-method family.

## License

Repository-level license: GNU Affero General Public License v3.0. See the repository root `LICENSE` and `LICENSE_NOTICE.md`.
