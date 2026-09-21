# Final scientific review of the astronomical-data addition

Review date: 2026-09-20. Scope: `paper/highz_analysis.tex`, `paper/highz_observed_table.tex`, relevant additions in `paper/main.tex`, and `reports/existing_data_results_and_outlook_cn.md`, checked against the current measurement, registration and design results. This review does not edit the manuscript or analysis scripts.

**Disposition: PASS for scientific interpretation as an exploratory existing-data analysis and observing outlook. No blocking scientific error remains.** The two requested notation/reporting clarifications have been implemented and independently re-read. The scientific interpretation remains conditional and does not turn the pilot into a cosmic-amplitude measurement.

## Direct consistency checks

| Manuscript/report statement | Check against present files | Finding |
|---|---|---|
| Representative equivalent widths: HE 1608, 2374, 2586 and Q0347 1608 | `observed_profile_metrics.csv`: 179.1980 ± 0.4990, 220.8385 ± 0.4356, 442.9815 ± 0.3776 and 235.5681 ± 0.6511 mÅ | Correct after stated rounding |
| HE centroid conditional statistical errors 18–117 m/s | Unrounded range 17.7689–117.0741 m/s | Correct |
| Prescribed continuum sensitivity 99–439 m/s | Unrounded range 99.0970–439.0320 m/s | Correct; text correctly distinguishes a prescribed sensitivity from measured systematic noise |
| Boundary sensitivity 0.81–1.66 km/s | Unrounded range 0.810843–1.660165 km/s | Correct; text correctly explains that changed intervals include different gas |
| HE 2382 and 2600 fractions below normalized flux 0.1 | 0.361905 and 0.317308 | Correct; treated as saturation warnings rather than precise optical-depth recovery |
| 2586 whole-window registration approximately +17 m/s and nominal residual ratio 5.26 | Expected-fluctuation/power configuration: +16.6601 m/s, 5.26291 | Correct; final prose identifies the error-array choice explicitly |
| 2586 red-window registration approximately −34 m/s and ratio 0.88 | Expected-fluctuation/power configuration: −34.2933 m/s, 0.880368 | Correct; not independent confirmation |
| 2600 whole-window residual ratio approximately 34.2 | Expected-fluctuation/power configuration: 34.1724 | Correct; withholding knot/boundary curvature errors is appropriate |
| Assumed logarithmic excess at the two epochs: 1.38% and 2.72% | `forecast_summary.json`: 0.0138351114 and 0.0272389907 | Correct for the declared background and Planck-time scale |
| Difference of normalized age basis approximately 0.495 | `data_informed_design.json`: 0.4945759366 | Correct; hypothetical B=100 m/s implies 49.4576 m/s mode difference |
| Exposure factors to native-pixel ratio 30: approximately 523 and 211 | Source coverage C/N=1.3118215 and 2.0633535; square-ratio factors 522.9896 and 211.3953 | Correct as conditional arithmetic, not telescope-time promises |

The continuum-to-error values in the old full coverage windows differ from those in the smaller new integration windows. The text correctly says **full-window** for the quoted 1.31 and 2.06 planning inputs, so these values should not be silently replaced by the smaller-window values in `observed_profile_metrics.csv`.

The moment equation agrees with the script's linear optical-velocity coordinate, signed flux deficits and partial-pixel weighting. It does not use the logarithmic velocity convention suggested as one possible choice in the earlier method notes. There is no inconsistency: the implemented linear convention makes the stated wavelength/velocity equivalent-width conversion exact within the piecewise-constant pixel representation.

## Interpretation and response audit

The manuscript makes the essential physical distinction between an epoch-dependent physical imprint and an earlier observer's computational ability. It does not claim that present-day reception of ancient light measures ancient observer resources. It also does not interpret six selected transitions as a complete Hamiltonian level sequence or as a GUE nearest-neighbor sample.

The new level-shift ansatz is explicitly described as **additional to the original standard-deviation schedule**. This resolves the main category risk in moving from noise precision to deterministic atomic shifts. The alternative line-noise response is appropriately said to require a temporal correlation function and line-shape model. No identification of the fitted empirical shift with alpha variation, cosmic noise amplitude or accessible zero height is made.

The rank example is mathematically consistent with its stated assumptions: six low-redshift line responses, one high-redshift line, free absorber velocities and a constant response-mode offset. The cosmic column adds no rank (four before and after); the analytic null vector reproduces this failure. The nuisance offset is an explicit assumption, so the prose correctly says the failure occurs **when that offset and absorber redshifts are free**. It should not be generalized to a claim that no single-epoch calibrated experiment can ever constrain a fully specified physical law.

The real-data section does not report a cosmic coefficient, cutoff limit, p=2 preference or a false null measurement. The distinction between an unidentified parameter and a measured zero is maintained in both languages. The main abstract and conclusion accurately summarize the new empirical contribution and its limits.

## Registration, conventional effects and source wording

The descriptive power transformation is correctly identified as insufficiently physical for precise transition-frequency inference. The noncommuting convolution relation, known Fe II 2586 additional absorption, unresolved saturation and conventional profile structure provide proper explanations for residuals without invoking new physics. The 0.88 subregion statistic is not treated as validating the entire model or as a calibrated goodness-of-fit probability. All 24 configurations are retained, including poor and boundary fits.

The discussion of omitted covariance is clear. A retained local curvature error remains conditional even if a fit avoids interpolation knots. The manuscript does not convert the large nominal residual statistics into significance against conventional atomic physics.

The published 2586 blend statement is supported by Murphy et al. (2022), sections 3.2 and 3.5, and the source files archived for this project. Saying that no numerical mask was fabricated is factually correct, although an optional more neutral manuscript formulation would be: “The exact published exclusion mask has not been transferred, so we retain the uncorrected profiles for this diagnostic.” This is an editorial suggestion, not a scientific blocker.

## Clarifications implemented and one optional citation refinement

1. **Implemented: named error array for the three headline registration values.** The English manuscript and Chinese report now identify the expected-fluctuation array for +17 m/s/5.26, −34 m/s/0.88 and the 2600 ratio 34.2. With the other source statistical-error array the corresponding residual ratios are approximately 6.08, 1.02 and 37.37. Both versions support the same interpretation.
2. **Implemented: complete transition-response equation.** The manuscript now gives `Delta nu_ul(t) = C (g_u − g_l) [h_2(t) − 1] / h_P` and defines h_P as Planck's constant. It also explicitly identifies the rank calculation's response vector as artificial and illustrative, rather than a derived Fe II response.
3. **Optional citation precision:** the six-line ESPRESSO redshift range is sound geometric arithmetic using the nominal 380–788 nm range from the existing instrument audit. An explicit citation to the ESO instrument-mode page would improve this paragraph's source traceability; Murphy et al. (2022) reports the actually used spectrum, whose wavelength range is not identical to the nominal range.

## Outlook audit

The design follows the observed constraints: weak and strong-line controls, common transitions across absorber redshifts, exposure-level calibration and modelling, sample selection frozen before confirmation, and independent instruments. It does not assume that downloading 467 spectra produces 467 eligible absorbers. It does not treat a quasar emission redshift as a usable absorber redshift.

The nominal ESPRESSO six-line range 1.36–2.03 is correct geometrically, with actual atmospheric/instrumental usability explicitly separate. The ANDES 0.4–1.8 μm, resolving power near 100,000 and 2035 first-light target are supported by the 2026 team paper and appropriately conditional. The 6–10 targets and 50–200 h are planning scenarios subject to target calculations and telescope availability. The manuscript does not promise discovery or an operational ANDES survey on a fixed date.

## Validation state at final check

During this review an older registration-validation snapshot was replaced following the root agent's global-search updates. The final file checked at **2026-09-20 13:54:17 UTC** reports **PASS**, contains no failed checks, and its `source_results_sha256` matches the current `observed_registration.json`. The previous failed snapshot is not the reviewed final result. That independent numerical validation is separate from the scientific limitations documented here.
