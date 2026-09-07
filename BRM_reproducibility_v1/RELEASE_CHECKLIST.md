# Release checklist for BRM v1.0.2

1. Confirm the final package metadata identify version 1.0.2 and the manuscript title exactly.
2. Verify the public OSF project URL used in the manuscript is https://osf.io/4u6dk/ and opens without authentication.
3. Verify the published GitHub Release is https://github.com/8iancachagasribeiro/_signal_pipeline-/releases/tag/brm-v1.0.2 and contains `BRM_reproducibility_v1_v1.0.2.zip`.
4. Confirm the release archive contains `mcphases_analyses.py`, `actigraphy_replication.py`, `ssf_estimators.py`, `empirical_audit/`, `make_empirical_figure.py`, and `validate_empirical_audit.py` in addition to the canonical simulation package.
5. Extract the release ZIP and run `python -m py_compile *.py`.
6. Run `python validate_brm_outputs.py --root .` and confirm `VALIDATION PASS`.
7. Run `python validate_empirical_audit.py --root .` and confirm `EMPIRICAL VALIDATION PASS`.
8. Run `sha256sum -c SHA256SUMS.txt` and confirm all packaged files pass.
9. Confirm `RELEASE_ASSET_SHA256.txt` matches the published release asset digest shown by GitHub.
10. Confirm the manuscript describes the ICC as `ICC(C,1)` consistency, not absolute agreement.
11. Confirm the manuscript and package describe internal interpolation only as an FFT-construction step and resampling only on originally observed days.
12. Confirm the manuscript cites only the v1.0.2 GitHub Release and public OSF project; no obsolete OSF, Zenodo, v1.0.0, or v1.0.1 submission destinations remain.
13. Optionally deposit the exact release ZIP in Zenodo. Add a Zenodo DOI to the manuscript only after the matching deposit is publicly accessible.
14. Re-run the manuscript link, metadata, numerical-anchor, accessibility, and full-page rendering audit before resubmission.