# Release and Zenodo checklist

1. Confirm the `brm-methodological-expansion` branch contains the final `BRM_reproducibility_v1/` code and metadata.
2. Download or use the locally validated `BRM_reproducibility_v1.zip` release asset.
3. Verify it with `sha256sum -c RELEASE_ASSET_SHA256.txt`.
4. Extract the ZIP and enter the extracted `BRM_reproducibility_v1/` directory.
5. Run `python -m py_compile *.py`.
6. Run `python validate_brm_outputs.py --root .` and confirm `VALIDATION PASS`.
7. Run `sha256sum -c SHA256SUMS.txt` and confirm every packaged file passes.
8. Create a GitHub release/tag named `brm-v1.0.0` from `brm-methodological-expansion`.
9. Attach `BRM_reproducibility_v1.zip` to that release.
10. Deposit that exact ZIP in Zenodo, verify metadata, and obtain the DOI.
11. Replace the manuscript repository placeholder with the permanent Zenodo DOI and optionally the GitHub release URL.
12. Re-run the BRM submission audit after inserting the DOI.
