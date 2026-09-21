# Independent ESPRESSO mask, covariance and gauge review

Reviewed 2026-09-20. This audit evaluates source decoding, selection, AR(1) algebra and local identifiability. It does **not** perform a fit or certify optimizer convergence. Numerical audit metadata and the exact source snapshot SHA256 are in `results/espresso_mask_covariance_review.json`. Calculations use the published parameter initialization and oversample = 3 to keep this independent check inexpensive.

## Source mask and region selection: PASS

The actual FITS primary array is 9 × 550806. Its header identifies row 0 as normalized flux, row 1 as normalized error, row 2 as normalized expected fluctuation, and row 4 as status, with 1 = valid and negative integers = clipped. The logarithmic wavelength convention uses a one-based reference pixel and is decoded correctly. The 0.4 km/s dispersion is represented internally by the corresponding logarithmic spacing.

The six source `fit_r_iso.f13` wavelength intervals contain 3033 pixels; 2931 remain after the stored status and positive-error mask. Counts, independently read from the actual FITS, are:

| Transition | Region pixels | Valid pixels | Stored exclusions |
|---|---:|---:|---|
| Fe II 2260 | 277 | 277 | None |
| Fe II 2344 | 515 | 430 | 22 pixels with status −5; 63 with status −7 |
| Fe II 2374 | 508 | 508 | None |
| Fe II 2382 | 579 | 579 | None |
| Fe II 2586 | 577 | 560 | 17 pixels with status −7 |
| Fe II 2600 | 577 | 577 | None |

Every selected line's error array is finite. The stored invalid pixels also have negative errors; none slips through with a positive error. These counts match the separately archived `feii_mask_audit.json`. The selection is an intersection of **author model region bounds** and **FITS-stored pixel exclusions**. They must not be described as the same object or as newly derived, data-blind blend masks. The −5 / −7 codes have not been assigned a physical cause by this audit. The author already inspected the data when choosing regions, components and exclusions; that limitation remains.

## Fixed AR(1) whitening: PASS

For each contiguous unmasked segment, the implemented transform retains the first standardized residual and sends later residuals to `(r_i − rho*r_(i−1))/sqrt(1−rho²)`. This is the stationary AR(1) whitening convention with correlation matrix `R_ij = rho^|i−j|`. It is applied consistently to every Jacobian column. The transform resets after every masked gap and at every transition.

Independent dense covariance matrices were built directly from the pixel indices, setting cross-segment covariance to zero. Solving the dense covariance equations gives agreement with the whitened quadratic form to relative error 0 at rho = 0.3 and 1.35e−16 at rho = 0.6. The largest absolute difference in the objective gradient contraction `J^T r` is 2.28e−12. This checks the actual selected data and actual Jacobian, rather than only a hand-generated example.

This is a correct **assumed covariance sensitivity model**. It does not estimate the actual resampling covariance or prove that the coadd noise is AR(1). Resetting after a masked gap also intentionally differs from retaining the very small stationary cross-gap correlations. For the observed minimum 17-pixel excluded run, the omitted rho = 0.3 nearest cross-gap correlation is only rho^18, but the segmentation convention must still be explicit.

## Gauge and conditioning: PASS with material inference limitation

With Fe II 2374 retained and its shift fixed, the five transition-shift directions remain locally identifiable after projection over the gas and continuum Jacobian. Their projected singular values at the published initialization are approximately 106.76, 106.06, 79.11, 26.75 and 6.80, in the code's velocity units. These numbers are numerical diagnostics, not final fitted confidence limits.

If the anchor were excluded while every remaining transition had a free shift, adding the same amount to all gas velocities and subtracting it from all line shifts would restore an exact gauge freedom. This issue was reported during review; the source now rejects free-shift models without the anchor. The guard was directly verified.

The full alternative Jacobian has 152 columns and rank 152 at the code's relative singular-value threshold 1e−10. Nevertheless, its condition number is approximately 7.35e8. **The original covariance implementation had a numerical inconsistency, now resolved.** It formed a pseudoinverse of `J.T @ J` with relative threshold 1e−12, approximately a threshold 1e−6 on J, and retained only 145 of the 152 directions at this initialization. It also squared the condition number by forming the normal matrix. These figures describe the historical implementation recorded in the original JSON audit metadata.

The revised source directly decomposes `J = U diag(s) Vh`, retains `s > s[0] * 1e−10`, and computes `covariance = (Vh[keep].T / s[keep]**2) @ Vh[keep]`. The cutoff now agrees with the reported Jacobian rank and the normal matrix is not formed. Direct source inspection confirms this correction; the follow-up source SHA256 is recorded in the JSON addendum. No new fit was run for this addendum.

The name `covariance_local_diagnostic` remains appropriate. A numerically consistent local covariance still does not justify ordinary Gaussian inference when nuisance directions are weak, parameters sit on bounds, or nonlinear alternatives have competing gas decompositions. Profiled likelihood or controlled refitting simulations are preferable for stronger uncertainty claims.

Optimizer success alone will not certify a unique gas decomposition or regular likelihood-ratio asymptotics. Final reports should retain optimum cost, gradient optimality, bounds and rank diagnostics, and describe nominal degrees of freedom as nominal. This audit has not inspected final fit output because no completed result was present when it ran.

## Disposition

Source masks and the AR(1) algebra pass. The missing-anchor gauge was guarded during review. The covariance threshold/normal-matrix issue has also been resolved by direct SVD with a consistent cutoff. Strong gas-parameter conditioning and bound-active local-Gaussian limitations remain and must constrain interpretation of any covariance or nominal significance. No mask or covariance implementation defect found here establishes an astrophysical signal, and no scientific null-hypothesis conclusion is claimed by this narrow audit.
