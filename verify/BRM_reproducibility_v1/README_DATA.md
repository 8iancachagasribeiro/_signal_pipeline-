# Data access and provenance

The Monte Carlo analyses in the archival release asset are fully reproducible without restricted human-participant data.

## mcPHASES

Raw mcPHASES data are not redistributed. The dataset is available through PhysioNet under credentialed-access conditions:

- Lin, B., Li, J. Y., Kalani, K., Truong, K. N., & Mariakakis, A. (2025). *mcPHASES: A dataset of physiological, hormonal, and self-reported events and symptoms for menstrual health tracking with wearables* (Version 1.0.0). PhysioNet. https://doi.org/10.13026/zx6a-2c81

Empirical reproduction begins from the original credentialed archive. The tagged repository source contains `mcphases_analyses.py` and `ssf_estimators.py`. The empirical pipeline treats each `study_interval` as a separate temporal segment, aggregates repeated eligible intervals back to one participant-level coupling estimate by precision-weighted Fisher transformation, and preserves the observed temporal calendar for surrogate analyses. Internal missing predictor days are linearly interpolated **only** to construct the regular grid required for Fourier phase randomization; surrogate values are subsequently evaluated only at originally observed days.

No identifiable participant data are included in the repository or release asset.

## Independent actigraphy datasets

The cross-domain actigraphy sources are publicly available from their original repositories:

- DEPRESJON: https://doi.org/10.5281/zenodo.1219550
- PSYKOSE: https://osf.io/dgjzu
- HYPERAKTIV: https://osf.io/3agwr

The tagged repository source contains `actigraphy_replication.py`. HYPERAKTIV diagnosis labels are read from source metadata (`patient_info.csv`); unresolved IDs are excluded rather than inferred from folder membership. FFT-based SSF calculations use regular, gap-aware temporal handling and do not silently compact missing rows.

The BRM materials do not redistribute those source datasets. Derived, non-identifying numerical audit summaries are retained where applicable for transparency and provenance.