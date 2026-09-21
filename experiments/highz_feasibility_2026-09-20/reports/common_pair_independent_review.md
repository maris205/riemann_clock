# Independent common Fe II pair uncertainty review

**PASS: 7,863 / 7,863 numerical checks; 160 new fit products; 8 actual optimizer-closure probes.** Both final refined summaries and their aggregate were audited. No optimization was performed by this review. All 160 scientific fit products report `ftol` termination (status 2); that is not a global-optimum certificate.

The observable is D = v(Fe II 2382) − v(Fe II 2374), in m/s. The 2374 shift is fixed at zero. A profile fit fixes only the 2382 shift; all remaining original H1 shifts and gas, continuum, zero-level and LSF parameters retain their original freedom and bounds. In particular, the D = 0 profile is a different model from H0, which fixes every relative line shift at zero.

## Final audited results

|Target|D (m/s)|Unrestricted local tangent sigma (m/s)|Profile objective rise 1 crossings (m/s)|Profile objective rise 3.84 crossings (m/s)|
|---|---:|---:|---|---|
|J232128-105122|-339.970553|179.581313|[-495.177706, -187.880853]|[-654.336850, -43.736228]|
|J004131-493611|-120.110048|399.412201|[-264.857331, 39.957201]|[-408.741608, 227.325730]|

These crossings are interpolated between actual optimized grid points. Levels 1 and 3.84 are descriptive objective contours, **not coverage-calibrated confidence intervals**. Active parameter bounds, selected gas architecture, diagonal expected-fluctuation weights, residual calibration and unknown pixel covariance limit their interpretation. In J004131, active narrow-component and instrumental-width boundaries make the profiled contour substantially tighter than the unrestricted tangent scale. Neither quantity should be advertised as physical precision.

|Target|Fairly restarted H0 chi2|H1 chi2|All-shift improvement|Added shifts|D = 0 one-coordinate objective rise|
|---|---:|---:|---:|---:|---:|
|J232128-105122|75.15888744|68.99884188|6.16004556|4|5.07438829|
|J004131-493611|95.22846063|87.25491207|7.97354856|2|0.58994264|

## Independent checks

- Recomputed the original endpoint local sigma using both a pivoted-QR nuisance projection and the full column-normalized Jacobian SVD covariance. For residual Jacobian J in km/s coordinates, sigma_D = 1000 / ||(I − Q Qᵀ) J_2382||. No reduced-chi-square scaling is applied. Original sigmas 179.582132 and 327.019043 m/s agree across independent methods to below 4e−13 m/s. Stable full ranks and singular-value-cutoff sensitivity are documented in `results/common_pair/review/projection_math.json`.
- Executed the actual optimizer closures with a recording substitute for the optimizer, for both targets at free D, D = 0 and D = ±1000 m/s. Verified exact parameter insertion, optimizer bounds, residuals and Jacobian column removal at initial and perturbed vectors. Temporary probe products were deleted; these were not scientific optimizations.
- Replayed every saved residual, full Jacobian, flux model, objective, parameter vector, data mask, native source index, expected-fluctuation error array, start mapping, bound-scaled gradient and active-bound list. Checked source/start hashes and every stopping flag. Baseline source and historical fit hashes remain unchanged.
- Independently checked finite differences for nine parameter classes at the original endpoints, and rechecked the actual final-basin shift, narrow-component log-width and LSF derivative. Verified source masks and row2 weighting against the already independently audited archive model.
- Checked endpoint inventories, minima across every recorded start, fixed-coordinate selection, final reference selection, local sigma, contour bracketing/interpolation, H0/H1 nesting, same-basin comparison arithmetic and exact aggregate/per-target summary equality. No lower recorded unfinished endpoint was ignored; all 160 recorded endpoints in this run terminated by ftol.

## Basin and metadata qualifications

J004131 first-pass profiling found a materially better nuisance basin: H1 chi2 changed from 104.53408492 to 87.25491207 and D from +25.635436 to −120.110048 m/s. The conventional H0 was restarted in that improved basin and reached 95.22846063. Comparing the old H0 to the new H1 would have been unfair; the exported comparison uses the new H0. The refined 21-coordinate J004131 curve must be used for interpolation. The union of historical and refined products spans 41 coordinates, including historical points trapped in the poor basin, and must not be plotted as one final optimized profile.

Both final references are fixed-coordinate profile endpoints, numerically below their best unrestricted fits by only 1.35e−9 and 1.05e−8 in chi2. They are equivalent within the explicitly recorded 1e−4 objective tolerance, but are not literally unrestricted optimizer outputs. The final point selection transparently retains this distinction.

The first-pass script completed and saved all 93 numerical fits, then its summary construction raised a duplicate `selection` metadata-key exception. A separate `common_pair_refine.py` reads those preserved outputs, performs the additional fits and successfully writes final summaries. This exporter exception did not corrupt the numerical products. The original script is preserved with its recorded hash and still has that documented summary-construction defect; this review does not claim that the first-pass script exits normally.

In `crossings`, the crossing itself is correct. However, `bracket_m_s` is sorted by coordinate while `bracket_delta_chi2` and `fit_files` are in inside-to-outside order. On the lower side these orders differ. Do not zip those arrays positionally; derive endpoint values from the referenced files or use named inside/outside endpoints.

The additional descending-order fits each start from the current best recorded basin, not from the preceding descending neighbor. Describe these as best-basin cross-starts. The contour routine returns the first crossing on each side of the selected minimum; it does not certify a complete, possibly disconnected likelihood set.

A future-run conditional weakness remains in the companion source: if the initial `refined_h1_from_h0` itself discovers a materially better basin but later profiles make no further improvement, the final-H0 restart condition may not trigger. This did not occur in the audited products. Future reuse should revisit that condition without silently changing hashes attached to the present experiment.

No result here establishes variation of alpha, a cosmic-time law, or physical uncertainty of mathematical zeros. Wavelength calibration, pixel covariance, architecture selection and broader model uncertainty are not marginalized. The closed contours and numerical convergence are conditional fitting diagnostics.

## Artifacts and reproducibility

- Main validator: `code/common_pair_validate.py` (runs no optimization).
- Numerical audit: `results/common_pair/review/validation.json`.
- Separate mathematical review: `results/common_pair/review/projection_math.json` and `.md`.
- Independent companion-source review: `results/common_pair/review/refinement_source_review.md`.
- Final scientific summaries: `results/common_pair/refined_summary.json` and each target’s `refined_summary.json`.

Run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_clock/experiments/highz_feasibility_2026-09-20/code/common_pair_validate.py
```

Final validator SHA-256: `e10733a9f8a8b3483e3444577759838823f30fa7b367194ba48b3b4d02ad583f`.
Final aggregate SHA-256: `fc91e39749a148a7a4b69d24e0f33895db3081e97abe8d9b4de926dd9059e857`.

## Final delivery review

**PASS: 126 / 126 additional delivery-only checks.** No new fits, model evaluations or whole numerical-audit rerun. Reviewed `common_pair_report.py`, the regenerated Chinese report, `harmonized_measurements.json`, export provenance and the final PNG visually. Fresh export/report/figure hashes are recorded in `results/common_pair/review/delivery_validation.json`.

The final figure uses only the selected summary profiles: 22 points for J232128 and 21 for J004131. It does not merge historical poor-basin points. Axes, target labels, D sign, tangent approximation, shading and the explicit non-calibrated-coverage caption are clear; there is no clipping of explanatory text. The chart intentionally focuses on objective rises up to 10 while the full profile remains in JSON.

All copied fields, source hashes, table rounding, H0/H1 comparisons, interval endpoints and bracket provenance match the audited final summaries. `D=v2382-v2374` is consistent with the kinematic log energy-ratio sign `-D/c`; that relation is not an alpha or cosmic-clock response. The schema explicitly excludes 2374 and 2382 from the list of other floating shifts, clarifies bracket ordering, and sets `precision_measurement=false` and `physical_time_inference=false`.

The Chinese report now qualifies the difference between the bounded nonlinear profile and unrestricted local approximation: hard boundaries may contribute, but nonlinear/residual-curvature effects have not been isolated. It correctly distinguishes D=0 profiling from the all-shift-zero H0, documents the newly found parameter basin and the first-stage exporter exception, and avoids calibrated-confidence or physical-time claims. No delivery blocker remains.
