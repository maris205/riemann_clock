# Independent numerical review of individual ESPRESSO exposure analysis

**PASS: 8,476 / 8,476 implementation and reproducibility checks; 0 failures.** The audit covers 459 saved line fits in nine configurations, 17 archived FITS files, the original UPL and atmospheric mask, and the paired diagnostic contrasts. This is not an audit of the truth of the conditional error model or the cosmic-time hypothesis.

The audited analysis source SHA256 is `87d26c6d5deb7d80f49b55efcf9eedbec47b834c25176643dd3804aa5d561f6f`. All final artifact hashes are recorded in `results/exposures/validation.json`; rerun `code/exposure_validate.py` after changing a scientific input or implementation.

## Independent reconstruction

- Verified each raw FITS file against both its manifest SHA256 and Git blob SHA1, and matched all 17 filenames to the original UPL exposure ordering.
- Reconstructed each selected pixel, vacuum-barycentric wavelength, midpoint boundary, flux, sigma, scale and velocity from the raw FITS. The independent nearest-coadd-pixel lookup uses full distance minimization; the fixed order masks are reparsed from the original UPL.
- Exactly reproduced all retained masks: 2,274 additional samples excluded by published manual clipping rectangles and 433 by the observer-frame atmospheric intervals transferred using exposure BERV and ±0.25 km/s padding.
- Midpoint bin widths agree with the supplied DLL widths to a maximum fractional difference of 1.4e-09 over the selected samples.
- Recomputed every saved model, residual and objective using QR variable projection instead of the fitting script’s least-squares solver. Independent Voigt-profile calculations reproduce the three intrinsic templates for both gas-parameter choices.
- Centered finite differences at two step sizes reproduce the local uncertainty matrices. The largest relative covariance-matrix difference from the saved result is 4.09e-05. An alternative fixed-beta nuisance-projected Fisher construction changes a shift standard error by at most 0.207%.
- Independently reconstructed the shared-anchor covariance, the weighted exposure means, heterogeneity statistics, and disjoint chronological contrasts. The nine early and eight late exposures are complete and non-overlapping.
- Verified that lower/upper order contrasts use identical 2374 anchor samples and estimates, so the common anchor cancels exactly. The two target-line order subsets have no common native rows. The trace comparisons have disjoint rows for all three lines. Independent paired-contrast propagation and aggregation reproduce the diagnostic files.

## Numerical refinement and limitations

Nine representative fits (exposure 0, three lines, primary/free-LSF/published-template configurations) were repeated with a 0.0125 km/s template grid instead of 0.025 km/s and 21 quadrature nodes instead of seven. The largest shift change was 5.03e-06 m/s, and the largest objective change was 5.41e-08. This checks numerical integration in the sampled cases; it is not a full optimization or model-uniqueness theorem.

All 459 saved fits have successful optimizer termination. 26 free-LSF fits touch a width bound (including repeated use of the same anchor in the lower/upper-order controls); none of the fixed-LSF primary fits touches a shift bound. Bound-active local covariance matrices are diagnostic curvature estimates, not validated Gaussian confidence intervals. The refinement and numerical replay do not remove this limitation.

The gas template is estimated from the same combined observations. Its uncertainty and shared calibration/extraction errors are not included in the reported covariance. Different photon samples are conditionally disjoint, but the shared template prevents treating their absolute means as independent astrophysical confirmation. The likelihood also assumes diagonal supplied ERRDATA; native sampling prevents final-coadd interpolation from being introduced here but does not establish intrinsically independent extraction or sky-subtraction errors.

The correct coadd comparison is the coadd fitted by the same fixed-template estimator. Comparing its result or the exposure result directly with the earlier full gas re-fit would confound the estimator and data-processing changes. The lower/upper-order discrepancy remains a calibration/extraction or model-response diagnostic, and the post-inspection p-value must not be presented as a calibrated new-physics discovery test.

Source-format details and review corrections are documented in `reports/exposure_provenance_review.md`.
