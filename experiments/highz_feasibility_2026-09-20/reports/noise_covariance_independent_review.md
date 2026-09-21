# Independent review of the ESPRESSO continuum-noise controls

Review date: 2026-09-21. Scope: `code/empirical_noise_controls.py`, saved continuum-control pixels and results in `results/noise_covariance`. This review does not repeat the nonlinear gas-model optimization or validate the forward model's physical completeness.

**Verdict: PASS for the saved numerical controls and their stated conditional interpretation. All 1,119 recorded numerical, provenance and scope checks passed. These are software/data consistency checks, not 1,119 independent scientific tests. No new discovery significance is established.**

## Independent reconstruction

The review reopened the source FITS file, reconstructed its wavelength grid from the header, independently parsed the isotope-weighted atomic wavelengths and all 35 fitting intervals in the author's three region files, and reconstructed all 36 candidate control tiles. It verified the first-owner assignment of source pixels, excluded regions, invalid pixels, source flux/error arrays, retained counts and all rejection reasons.

All saved control pixels are unique across tiles. The primary selection retains 25 disjoint tiles and 13,750 unique source pixels; none falls inside the excluded author fitting intervals. First-owner assignment occurs before validity and region exclusions, as documented. Consequently the geometric sensitivity cannot be interpreted as a search over all possible alternate pixel ownership schemes.

The linear continuum fits were recomputed using weighted `numpy.polyfit`, independently of the production polynomial-matrix solver. Lag products were independently calculated through a source-index lookup: both endpoints and every intervening source pixel must be retained. This reproduces the production rule `idx[i+k]-idx[i] == k`, avoids artificial adjacency across gaps, and correctly counts each within-tile lag pair once. Constant and quadratic detrending moments were also checked.

The per-tile moments, aggregate covariance/ACF, all transition and side splits, threshold variants, geometric-retention variants, median-flux-gate variants, and the 10,000 fixed-seed interval-bootstrap draws reproduce the saved results within the recorded numerical tolerances. All recorded code/data/baseline dependency hashes match the files used for this review.

## Reproduced continuum results

| Quantity | Independently verified result |
|---|---:|
| Selected tiles / pixels | 25 / 13,750 |
| Standardized residual variance | 0.7047517138 |
| Standardized residual RMS | 0.8394949159 |
| Lag-1 ACF | 0.3726390014 |
| Lag-2 ACF | 0.0440006825 |
| Lag-3 ACF | 0.0319338789 |
| Conditional cluster-bootstrap lag-1 interval | [0.3533045307, 0.3915642485] |

The ACF is a pair-averaged residual product normalized by the aggregate zero-lag residual variance. It is not a sample Pearson correlation separately centered for each lag, and the variance has no degrees-of-freedom correction after fitting two continuum parameters per tile. The output and report now disclose this distinction.

All 12 broad-bin/RMS gate variants retain the same 25 tiles. That fact alone is not broad selection robustness: it does not vary tile placement, ownership or geometric coverage. The additional retention-fraction scan includes 25–28 tiles and gives lag-1 ACF 0.37264–0.37575. The strict [0.99, 1.01] median-flux gate retains 24 tiles with lag-1 ACF 0.37587. These are useful bounded sensitivity checks, and do not eliminate post-selection uncertainty.

## Covariance and local-shift projection

One baseline forward-model/Jacobian evaluation was reproduced without running an optimizer. The fixed baseline contains 2,931 fitted pixels. The reviewer independently reconstructed nine contiguous valid-pixel segments across the six transitions. Covariance resets at every masked gap and at each transition boundary; no correlation is introduced between these segments.

For AR(1), an analytic innovation transform was used instead of production Toeplitz-Cholesky whitening. For MA(1), diagonal and Bartlett-tapered ACF cases, covariance matrices were independently assembled from lag distances and Cholesky-whitened. All matrices are positive definite for the reported cases. The covariance minimum eigenvalues and weighted residual norms reproduce the saved output.

The nuisance basis was independently calculated using SciPy's GESVD driver; five shift estimates were then solved using GELSD least squares. This differs from the production NumPy SVD projection and inverse implementation. All ten covariance cases at each of three singular-value cutoffs reproduce the stored ranks, estimates, covariance matrices and local score values. Each case omits a chi-square-tail probability in the final saved sensitivity output.

At relative cutoff 1e-10, the independently reproduced **local tangent-space scores**, not nonlinear likelihood-ratio improvements, are:

| Assumed covariance | Nominal error scale | Continuum-derived variance scale |
|---|---:|---:|
| Diagonal | 58.80782 | — |
| AR(1), measured lag-1 parameter | 28.09642 | 39.86711 |
| MA(1), measured lag-1 parameter | 33.93384 | 48.15007 |
| Tapered measured lag-0–10 ACF | 32.74208 | 46.45902 |

The uniform variance rescaling behaves exactly as it should: estimates remain invariant, their covariance multiplies by 0.7047517, and the local score divides by that factor. This independently confirms that the apparent partial recovery of score in the scaled cases is the expected algebraic consequence of shrinking assumed errors, not a separate signal detection.

## Scientific limitations that must remain attached

- The tiles and thresholds were defined after previous inspection. The results are exploratory; the bootstrap does not include selection uncertainty or testing multiple possible control placements.
- Known fitting intervals are excluded, but weak unrelated absorption, residual sky features, extraction artifacts and continuum structure can remain. A measured continuum ACF establishes correlated residual structure, without uniquely identifying its cause.
- The control noise is not demonstrated to be stationary across wavelength. Transition-group variances range from 0.67330 to 0.79583 and lag-1 ACFs from 0.34357 to 0.39996; only two selected tiles surround 2600. The common aggregate is a conditional choice.
- Disjoint wavelength tiles are not proven independent. Some are adjacent, and common calibration/continuum errors may span several tiles. The cluster-bootstrap percentile interval is therefore conditional and does not establish a complete uncertainty budget.
- Continuum normalized-error behavior does not automatically transfer into saturated absorption. The `control_scaled` rows must remain sensitivities, not preferred corrected error bars. A photon/background/exposure-level covariance analysis would be required for that inference.
- Lag-2 ACF differs substantially from the AR(1) extrapolation of lag-1 ACF. AR(1), MA(1) and tapered ACF are alternative approximations; none is established as the instrumental covariance process. Bartlett tapering changes the applied correlation relative to the raw measurement.
- The projected nuisance directions are unconstrained and can include weak or bound-active parameter combinations. The local score does not equal an achievable bounded nonlinear improvement and cannot replace nonlinear covariance-aware refits, gas-model selection or calibration controls. The report states this explicitly.
- The surviving conditional relative shifts neither identify varying constants nor establish a cosmic-time dependence. This one absorber cannot select an inverse-log-squared time law.

Machine-readable checks and the exact input hashes are saved in `results/noise_covariance/independent_validation.json`. The temporary independent audit implementation used for this review is `/tmp/audit_clock_noise.py`; the project artifact is the recorded checks/results, while the production experiment remains reproducible through its published script.
