# Final scientific-scope and delivery review: ESPRESSO conventional baseline

Date: 2026-09-20. Reviewer task: independent source/method and manuscript review. No new fits were run for this review; the PDF layout is handled separately by the main task.

**Verdict: PASS.** The minor derived-metadata correction is complete, and the baseline-bound clarification has been added to the manuscript. Neither changes fitted numerical results.

## Material inspected

- `paper/espresso_null_section.tex`, the generated numbers and comparison table;
- `paper/main.tex` abstract, motivation, limitations, conclusion and data declaration;
- `reports/conventional_null_results_cn.md`;
- the English and Chinese project README additions and experiment README replay instructions;
- selected-pair comparison CSV, local-control and velocity-split reports;
- initialization/continuation provenance and the saved numerical-validation summaries.

## Result and conditioning checks

The manuscript and reports consistently use the selected nonlinear comparison, H0 `null_cross.json` versus H1 `alternative_from_published.json`: chi²=2305.49054095 versus2253.77730148, improvement51.71323947 for five extra shifts. The baseline has2931 valid pixels and147 nominal parameters, giving2784 nominal degrees of freedom. The expected-fluctuation, AR(1), LSF, blend-exclusion and exploratory controls agree with the machine-readable selected comparison table. All18 selected fits across nine pairs have successful stopping flags; unselected budget-limited attempts remain documented. This is correctly distinguished from a certificate of global optimality.

The strong-line shifts, roughly−118 and−109m/s relative to2374, are consistently described as parameters of a restricted absorption fit. Their persistence under the implemented controls is reported honestly: the manuscript does not say that those controls explain the residual away. Removing2382/2600 yields improvement7.08993 for three shifts; reduced information is explicitly acknowledged, and this deletion is not presented as identification of a physical mechanism.

The local score58.80782 is kept separate from the nonlinear improvement51.71324. The20,000 simulations operate in a fixed Gaussian tangent space, not20,000 nonlinear nuisance refits. Exact-forward injection followed by local recovery is described accurately. Active bounds, weak nuisance directions, local steps outside the allowed domain and unknown pixel covariance are disclosed. Conditional local tail probabilities remain in the detailed technical report, with explicit language preventing their interpretation as new-physics significance.

The velocity split is correctly called post-inspection and local; compatible blue/red shifts are not promoted into independent replication. The two portions share one absorber and observational/calibration systematics.

## Scientific interpretation and sources

The source is the pinned historical ESPRESSO coadd of the same low-redshift absorber, not a2026 exposure or an additional cosmological epoch. The six FeII transitions differ from the earlier UVES set. The source paper's published free-alpha result is not treated as an exact fixed-alpha null. Both hypotheses refit the shared gas and continuum parameters; opacity, isotope weighting, intensity convolution, pixel averaging and zero-level corrections are stated consistently with the implemented model.

The description does not assume that any residual must be beyond standard astrophysics. Alternative component structures, saturation/covering, unrecognized blends, atomic data/isotopes, exposure-level calibration, instrumental kernel and full covariance remain unresolved explanations. The paper infers neither alpha variation nor a temporal law nor a Riemann-zero cutoff. This is the appropriate answer to the user's first-stage request: a real restricted physical baseline comparison has now been performed, and it has not completed the stronger task of excluding all conventional explanations.

The source paper'sFigure10 point estimates qualitatively agree with the present shift pattern after changing the reference fromMgI2852 toFeII2374. The project acknowledges that reuse of the same observations is not an independent discovery. The earlier source audit's visual numbers are explicitly approximate and are not inputs to the new test.

The64 forward-model checks and550 saved-point checks are reported as numerical verification, not physical validation or external peer review. The underlying validation summary files both reportPASS. I did not repeat their calculations in this final reading review.

## Resolved minor corrections/clarifications

1. The auxiliary derived file `data/raw/espresso_null/fit_r_iso_final_parsed.json`, created during the source audit, originally mislabeled the two numeric fields followingalpha as `temperature` and `coverage`. In the VPFIT `novars4` format these are the fixed turbulent contribution and fixed temperature; the1.00E+00 input is1K. The source rawf18 is unchanged and the fitted code reads only its N,z,b columns directly, so this metadata error has no effect on any fit. The derived labels have now been corrected to `bturb_km_s` and `temperature_K`, with matching tag names, for all191 ionic rows. The correction is recorded in the method audit; no result rerun was needed.
2. The earlier manuscript mentioned bounds without spelling out the baseline ranges: velocities within±2km/s of the published initialization, b=0.15–30km/s, relative shifts within±0.5km/s, continuum coefficients within±0.03, and column densities within±2dex subject to absolute logN bounds7–17. These ranges have now been added to the manuscript. The additional observed adjacent-residual correlation is explicitly identified as a diagnostic of residuals, not a measurement of pure-noise covariance.

No other correction is required to the scientific conclusion reviewed here. This review is not a recommendation to launch further fitting before delivering the completed, carefully scoped result.
