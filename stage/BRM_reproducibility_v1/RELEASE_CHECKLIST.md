# Release checklist for BRM v1.0.1

1. Confirm the published GitHub Release is https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.1 and is neither a draft nor a prerelease.
2. Confirm the validated archival asset is named `BRM_reproducibility_v1_v1.0.1.zip`.
3. Confirm the release asset SHA-256 exactly matches `RELEASE_ASSET_SHA256.txt` and the checksum shown in the release description.
4. Use the **release page and release tag**, not a mutable branch-tree URL, as the stable repository citation target.
5. Verify the public OSF project URL used in the manuscript is https://osf.io/4u6dk/ and opens without authentication.
6. Extract the release ZIP and run `python -m py_compile *.py`.
7. Run `python validate_brm_outputs.py --root .` and confirm `VALIDATION PASS`.
8. Run the package checksum validation and confirm all packaged files pass.
9. Confirm the tagged repository source contains the audited empirical scripts `mcphases_analyses.py`, `actigraphy_replication.py`, and `ssf_estimators.py`; raw mcPHASES data must remain absent.
10. Confirm manuscript terminology matches the implemented statistics: recovery fidelity for ordering, **ICC(C,1) consistency ICC** for cross-person consistency, and RMSE/bias for magnitude error.
11. Confirm the manuscript contains no obsolete OSF/Zenodo destinations and makes no formal preregistration claim.
12. Optionally deposit the exact release ZIP in Zenodo. Add a Zenodo DOI only after the deposit is published and independently accessible.
13. Re-run manuscript link, metadata, numerical-anchor, accessibility, and visual-rendering audits immediately before resubmission.