# Corrected release checklist for BRM v1.0.1

1. Confirm `brm-v1.0.1-candidate` contains the final `BRM_reproducibility_v1/` code and corrected metadata.
2. Verify the OSF project URL used in the manuscript is https://osf.io/4u6dk/ and is accessible to the intended audience.
3. Do not describe the OSF project as a formal preregistration unless a separate timestamped, immutable OSF Registration is documented.
4. Build a fresh `BRM_reproducibility_v1.zip` from the corrected candidate contents.
5. Generate a new SHA-256 for that exact ZIP and update `RELEASE_ASSET_SHA256.txt` if the archive contents change.
6. Extract the ZIP and run `python -m py_compile *.py`.
7. Run `python validate_brm_outputs.py --root .` and confirm `VALIDATION PASS`.
8. Run the package checksum validation and confirm all files pass.
9. Create a GitHub release/tag named `brm-v1.0.1` from the corrected candidate commit.
10. Attach the freshly validated ZIP to that release.
11. Optionally deposit that exact ZIP in Zenodo. Add a Zenodo DOI to the manuscript only after the deposit is published and the DOI resolves publicly without authentication.
12. Update manuscript repository/version references to the final `brm-v1.0.1` release URL.
13. Re-run the manuscript link, metadata, numerical-anchor, and rendering audit before resubmission.