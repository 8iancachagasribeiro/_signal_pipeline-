# Changelog

## 2026-08-24 — methodological audit branch

### Blocking corrections
- Removed the use of empirical smooth-signal fraction (SSF) as if it were classical reliability in power simulations.
- Added `ssf_power.py` to calibrate simulated predictor and outcome series directly to target spectral SSF.
- Added explicit `coupled_fraction` (`q`) because SSF does not identify how much smooth outcome variance is E2-coupled.
- Reframed registered-test power and budget analyses as sensitivity surfaces rather than point estimates of empirical power.
- Added participant-specific-predictor random-slope LRT in `fastlrt.py`.
- Changed mcPHASES surrogate generation to preserve the calendar before Fourier phase randomization.
- Changed group summaries in the differential-prediction analysis to person-centered within-person associations.
- Replaced direct SSF-based disattenuation claims with an explicit `f` sensitivity analysis.

### Reproducibility corrections
- Submission-grade theoretical grids default to 500 Monte Carlo replicates where applicable.
- Removed environment-specific output paths.
- Added portable result paths across scripts.
- Corrected HYPERAKTIV labeling to use `patient_info.csv`; unresolved labels are never guessed.
- Actigraphy analyses now save tables used to determine current sample sizes.
- Rebuilt `make_figures.py` to implement all five current Article 1 figures with current numbering.
- Added `audit_consistency.py` for static and numerical smoke testing.
- Updated README so claim-to-code mappings no longer cite superseded power numbers or stale figure/table numbering.

### Consequence for the manuscript
The masking mechanism remains reproducible. Exact claims that empirical SSF implies a fixed attenuation, zero power, or a specific gain from repairing instruments must be replaced by the audited sensitivity analyses before submission.
