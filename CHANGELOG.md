# Changelog

## 2026-08-24 - full methodological and empirical audit

### Blocking methodological corrections
- Removed the use of empirical smooth-signal fraction (SSF) as if it were classical reliability in power simulations.
- Added `ssf_power.py` to calibrate simulated predictor and outcome series directly to target spectral SSF.
- Added explicit `coupled_fraction` (`q`) because SSF does not identify how much smooth outcome variance is E2-coupled.
- Reframed registered-test power and budget analyses as sensitivity surfaces rather than point estimates of empirical power.
- Added participant-specific-predictor random-slope LRT in `fastlrt.py`.
- Replaced direct SSF-based disattenuation claims with an explicit identification-sensitivity analysis.

### mcPHASES empirical re-audit from the original 1.0.0 archive
- The dataset contains 42 unique participants and 62 participant-intervals. Each `study_interval` is now a separate temporal segment for FFT and surrogate operations; the long 2022-to-2024 gap is never bridged.
- Repeated intervals are Fisher-z aggregated back to one participant-level coupling estimate, so intervals are not treated as independent participants.
- Pooled within-person summaries are centered within participant-interval.
- Fitbit files described as daily are collapsed to one participant-interval-day before joining to hormone data. Resting-heart-rate values <=0 are treated as missing.
- Confirmatory heterogeneity remains null after correction: fatigue p=.5329; mood swing p=.8802.
- The previously reported exploratory skin-temperature result was reproduced as a legacy raw duplicate-merge artifact (p approximately .024) and disappears after daily de-duplication plus episode-aware analysis (p=.8743). Resting heart rate is also null (p=.1557).
- Corrected participant-balanced spectral SSF estimates are E3G=.4271, fatigue=.3618, mood swing=.3783, cramps=.5724, and bloating=.3491. Objective SSF now travels with eligible-participant counts because strict regular-grid gap handling changes the estimand materially.
- Cramping remains the strongest phase-locked result after within-interval standardization: |r| with E3G level=.0782, phase eta-squared=.1744, menstrual mean=+0.836 within-interval SD.
- Added `audit_results/2026-08-24/` containing derived, non-identifying outputs from the credentialed source archive.

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
The masking mechanism remains reproducible. The manuscript must remove the old skin-temperature p approximately .027, the claim that objective temperature reveals heterogeneous E3G coupling where self-report does not, and any statement that SSF directly identifies reliability, empirical power, or exact gains from instrument repair. The empirical mcPHASES section should instead report the corrected null surrogate results, participant-interval handling, revised SSF estimates, and the retained phase-misspecification result for cramping.
