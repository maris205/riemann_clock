## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate, bounded independent fit review
- Origin Date: 2026-09-21
- Verification Status: VERIFIED for the scoped input, objective and selection checks; scientific interpretation remains conditional.
- Version Label: j0530_fit_review_v1

# J053007−250329: independent review of the completed fitting campaign

The saved results support a modest conditional improvement when three relative line shifts are added. They do not establish anomalous atomic physics, a change in α, or a cosmic-time law. The updated Chinese result report and the actual profile figure present this distinction adequately. No changes to fit code, fitted values, selection, or masks were made during this review.

## Numerical reproduction

**353/353 scoped checks passed.** The companion file is `results/archive_complete/J053007-250329/independent_fit_review.json`.

The selected H0 and H1, the H0 cross-start, and all five fixed-shift profile fits use exactly the same **480 native pixels**, comprising 120 pixels each for Fe II 1608, 1611, 2374 and 2382. Their source-index lists match the independently reconstructed pre-fit selection. Fluxes, wavelength coordinates, and expected-fluctuation errors reproduce the full coadd without a new mask. All parameters lie within the recorded bounds. H1 adds only the three recorded relative shifts to the same nuisance-parameter set.

For these eight saved fits, I independently reconstructed each normalized residual from its saved physical profile and the original full-coadd flux/error arrays. The residuals, total and per-line objectives, nominal parameter counts, active-bound reports, and AICc arithmetic reproduce the saved values. Re-evaluating the current shared model at the saved parameters reproduces the objectives and Jacobians. This is a reproducibility replay, not a second implementation of the complete Voigt/convolution calculation; the separate general numerical audit remains relevant.

| Quantity | Reproduced result |
|---|---:|
| H0 χ² / nominal degrees of freedom | 482.2984656 / 416 |
| H1 χ² / nominal degrees of freedom | 478.1224282 / 413 |
| Additional shift parameters | 3 |
| H0 minus H1 χ² | 4.1760374 |
| H0 AICc | 630.3466583 |
| H1 AICc | 634.2389330 |
| H1 minus H0 AICc | +3.8922747 |

Thus the smaller raw objective under H1 does not translate into a preference for H1 under this same conditional AICc bookkeeping. AICc itself remains conditional on the diagonal noise model, selected architecture, and nominal parameter accounting; it is not a general calibration of evidence in the presence of active bounds or gas degeneracies.

## Architecture and optimization

The recorded component counts are 4, 6, 8, 10, 12 and 16. For each count, the summary selects the best saved successfully stopped H0 among its finite initializations. The minimum H0 AICc occurs at 16 components; no free-shift objective enters that component selection. The implementation completes this H0 selection and its operational adequacy check before running H1. This verifies the recorded execution design, not an independently certified chronology of every possible external analysis.

The selected count is the **largest tested count**. AICc has not reached a demonstrated plateau over a broader grid, and a unique gas decomposition has not been established. The sharp improvement from insufficient component sets to the 16-component set is a useful warning against interpreting residuals from sparse gas models as exotic frequency shifts.

Both selected fits stopped by **ftol**, with reproduced first-order optimalities approximately **0.200** for H0 and **0.640** for H1. Those are not values satisfying the requested gradient tolerance of 10⁻⁶. The successful stop flag means the documented stopping rule fired; it does not prove a stationary point or global optimum. The close H0 cross-start objective and sparse-profile center are useful local stability checks but do not exhaust alternative solutions.

Both fits place `zero_1611` and `zero_2374` at their **+0.02 upper limits**. No selected relative-shift or LSF parameter is reported at its bound. Active zero-level limits are especially relevant because the strong-line cores are saturated and the weak 1611 line has less independent leverage on its zero level. The quoted local error does not propagate changing those bounds or changing the component architecture.

## Shift curvature and actual profiles

The saved Fe II 2382 minus 2374 shift is **+268.781 m/s**, with local conditional standard deviation **130.094 m/s**. An independent pivoted-QR nuisance projection reproduces the stored SVD result and nuisance rank 66. The projection includes every other gas, continuum, zero, LSF, and shift column. It remains a local linear result that ignores active-bound truncation, calibration uncertainty, blends, architecture uncertainty and unknown pixel covariance.

The five nuisance-reoptimized profile points have objective increments approximately 3.907, 0.976, −0.000008, 0.974 and 3.890 at the local center plus −2, −1, 0, +1 and +2 nominal standard deviations. This resembles the expected local quadratic shape; it does not establish interval coverage. The central point improves the free objective by only **8.03 × 10⁻⁶**, which is retained rather than silently treated as zero. No materially lower profiled solution appears on this sparse grid.

I inspected `physical_profiles.png`. The figure displays all four actual flux profiles, both model curves, error bars, and separate normalized-residual panels. H0 and H1 are nearly superposed, consistent with the limited objective improvement. The 2374 line still has visibly structured departures, particularly on the blue shoulder and portions of the red wing; its H0 χ² per pixel is **1.625**, larger than those of the other three lines. Passing the recorded global and per-line thresholds is an operational permission for a conditional comparison, not evidence that the residuals are white or that all conventional effects are absent.

The separate quadrature file reports that fixed-parameter refinement from 21 to 55 integration subdivisions changes the comparison by **1.43 × 10⁻⁶**. Its OS21 objectives match the selected pair. That report assesses integration at fixed parameters; it does not replace optimization or gas-architecture checks.

## Scientific wording and fallacy scan

The result report now explicitly states the finite initialization search, highest tested component count, non-stationary stopping limitation, active zero bounds, saturation, and unpropagated covariance/calibration uncertainties. These qualifications address the material issues raised during review. It is suitable as an exploratory conventional-model result, provided later summaries preserve those limits and do not call the archive sample exhaustive or the fit globally converged.

11/11 statistical-fallacy categories were considered. Relevant concerns are sample selection/survivorship, look-elsewhere effects, and researcher choices in windows, transitions, component sets and noise models. The fixed current comparison and saved unsuccessful endpoints improve transparency but do not turn an exploratory analysis into a confirmatory detection. Group-level causal claims, ecological inference, collider adjustment, prevalence inference, regression-to-the-mean inference, reverse causality and Simpson aggregation claims are not made here. In particular, a conditional displacement is not causally attributed to cosmological evolution.

**Review disposition: PASS WITH LIMITATIONS.** The arithmetic, pixel identity, conditional selection, and reported local curvature reproduce. The scientific conclusion should remain that this absorber does not add compelling evidence for extra relative shifts under the tested model family.

## Final report verification

The finalized Chinese report adds the reproduced AICc comparison and links to the 38-check input review and 353-check fit review. These additions are accurate and retain the conditional interpretation. Eight finalization checks passed: the AICc values/difference, review references/counts, unchanged figure PNG hash, and unchanged hashes of all eight reviewed fit JSON/NPZ pairs. The reviewed report hash is now `90ff9695215642d559738209bf50f06c717c33cfd518f923e2c24337e5723d35`. No fit or mask changed, and the original 353-check count is retained.
