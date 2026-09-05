# Release checklist for BRM v1.0.1

1. Confirm `brm-v1.0.1` contains the final `BRM_reproducibility_v1/` code and corrected metadata.
2. Verify the public OSF project URL used in the manuscript is https://osf.io/4u6dk/ and opens without authentication.
3. Verify the published GitHub Release is https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.1 and contains the validated archival ZIP.
4. Verify the canonical source snapshot is https://github.com/8iancachagasribeiro/_signal_pipeline-/tree/brm-v1.0.1/BRM_reproducibility_v1.
5. Extract the release ZIP and run `python -m py_compile *.py`.
6. Run `python validate_brm_outputs.py --root .` and confirm `VALIDATION PASS`.
7. Run the package checksum validation and confirm all packaged files pass.
8. Confirm `RELEASE_ASSET_SHA256.txt` matches the currently published release asset.
9. Optionally deposit that exact ZIP in Zenodo. Add a Zenodo DOI to the manuscript only after the deposit is published and the DOI resolves publicly without authentication.
10. Re-run the manuscript link, metadata, numerical-anchor, and rendering audit before resubmission.