# Canonical full Monte Carlo outputs

The compact summary tables required to inspect the principal BRM numerical anchors are kept in Git. The full canonical CSV outputs and publication-ready figure files are bundled in the published v1.0.1 archival release asset:

https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.1

The corresponding source snapshot is:

https://github.com/8iancachagasribeiro/_signal_pipeline-/tree/brm-v1.0.1/BRM_reproducibility_v1

After extraction, validate with:

```bash
python validate_brm_outputs.py --root .
sha256sum -c SHA256SUMS.txt
```

The published release archive contains:

- `results/brm_robustness_metrics.csv` (400 cells)
- `results/brm_sampling_design.csv` (384 cells)
- `results/brm_ssf_benchmark.csv` (576 method rows)
- `results/brm_surrogate_calibration_power.csv` (31 cells)
- six publication-ready SVG figures
- executable scripts, release metadata, seeds, environment record, and compact tables.

The current outer archive checksum is recorded in `../RELEASE_ASSET_SHA256.txt`.