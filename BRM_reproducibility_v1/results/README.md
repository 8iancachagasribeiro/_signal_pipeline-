# Canonical full Monte Carlo outputs

The full canonical CSV outputs are distributed in the validated release asset `BRM_reproducibility_v1.zip` rather than duplicated in Git history.

Verify the release asset against `../RELEASE_ASSET_SHA256.txt`, extract it, and run:

```bash
python validate_brm_outputs.py --root .
sha256sum -c SHA256SUMS.txt
```

The release archive contains:

- `results/brm_robustness_metrics.csv` (400 cells)
- `results/brm_sampling_design.csv` (384 cells)
- `results/brm_ssf_benchmark.csv` (576 method rows)
- `results/brm_surrogate_calibration_power.csv` (31 cells)
- six publication-ready SVG figures
- executable scripts, release metadata, seeds, environment record, and compact tables.
