# Release checklist for BRM v1.0.1

1. Confirm `brm-v1.0.1` contains the final `BRM_reproducibility_v1/` code and corrected metadata.
2. Verify the public OSF analysis-plan URL used in the manuscript is https://osf.io/4u6dk/ and opens without authentication.
3. Do not describe the OSF project as a formal preregistration unless a separate timestamped, immutable OSF Registration is documented.
4. Use the canonical repository snapshot https://github.com/8iancachagasribeiro/_signal_pipeline-/tree/brm-v1.0.1/BRM_reproducibility_v1 in submission materials until a matching GitHub Release or DOI archive is published.
5. Build a fresh `BRM_reproducibility_v1.zip` from the `brm-v1.0.1` contents before any archival release.
6. Generate a new SHA-256 for that exact ZIP and replace the pending marker in `RELEASE_ASSET_SHA256.txt`.
7. Extract the ZIP and run `python -m py_compile *.py`.
8. Run `python validate_brm_outputs.py --root .` and confirm `VALIDATION PASS`.
9. Run the package checksum validation and confirm all packaged files pass.
10. If creating a GitHub Release, use tag `brm-v1.0.1` and attach the freshly validated ZIP.
11. Optionally deposit that exact ZIP in Zenodo. Add a Zenodo DOI to the manuscript only after the deposit is published and the DOI resolves publicly without authentication.
12. If a GitHub Release is published, update manuscript repository references from the versioned branch snapshot to https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.1.
13. Re-run the manuscript link, metadata, numerical-anchor, and rendering audit before resubmission.