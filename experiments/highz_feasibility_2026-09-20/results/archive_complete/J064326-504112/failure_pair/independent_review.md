# Independent review of the J0643 failure-case shift comparison

**93 / 93 grouped checks passed.** The newly authorized comparison is reproducible and correctly remains a failed-model diagnostic. No optimizer was run during this review. Checks compare scientific quantities or whole array groups; individual pixels are not counted as separate tests.

The original quality gate was a sample-quality screen, not a test establishing population-level constancy. An inadequate ordinary model could itself contain sensitivity to a physical frequency change. It is therefore appropriate to report the explicitly amended, bounded H1 comparison, including its failure to repair the residual structure, rather than infer a null population result from the original exclusion.

## Sequence and preservation

Exactly two new results exist: one three-shift H1 initialized from the protected 16-component primary and one H0 cross-start initialized from that H1. Both starting vectors match the prescribed parameter mapping and clipping. The three added shifts correspond to Fe II 1608, 1611 and 2382, with 2374 anchored at zero and all shifts bounded at ±1 km/s. The H0 model is reproduced exactly inside H1 when those added shifts are zero.

The H0 cross-start gives χ² = 532.225296838, slightly worse than the original primary's 532.225273754. The 0.01 improvement trigger for another H1 was not met. Retaining the original H0 is correct, and the summary's selected-null path resolves to the protected primary record outside `failure_pair`.

All 40 original fit files—20 JSON records and their 20 NPZ arrays—remain hash-identical to the earlier independent audit. The original primary, configuration, runner and engine hashes match the amendment. Both new fits use the same current coadd, metadata, atomic table, base-model and adapter hashes. The amendment timestamp precedes the saved fit results, and inspection of the runner confirms that it writes the amendment before calling the fit engine.

All four lines retain the exact same 436 native pixels, velocities, wavelengths, fluxes, statistical errors, expected-fluctuation weights, masks and source indices as the primary. Direct reconstruction from the full coadd confirms the fixed window and valid-pixel selection. No pixel rejection or gas-component reselection occurred.

## Recomputed comparison

| Quantity | H0 | H1 |
|---|---:|---:|
| Gas components | 16 | 16 |
| Native pixels | 436 | 436 |
| Fitted parameters | 64 | 67 |
| χ² | 532.225273754 | 527.587722694 |
| Nominal degrees of freedom | 372 | 369 |
| Fe II 1611 χ² per pixel | 2.058481575 | 2.031517966 |
| Original conditional residual gate | Fails | Fails |

The independently recomputed improvement is **Δχ² = 4.637551060 for three added parameters**. Both models pass the aggregate χ²/ν threshold of 1.5, but 1611 alone exceeds the per-line χ²/pixel threshold of 1.8. The bounded shifts therefore do not repair the principal residual inadequacy in this model.

H1's conditional fit coordinates are +5.55 m/s for 1608, −895.27 m/s for 1611 and −217.02 m/s for 2382. These are diagnostic coordinates of an inadequate model, not calibrated displacement measurements. The 1611 shift is near the lower permitted range but does not activate its bound. Both selected models retain the `logb_14` and `zero_1611` active bounds. H1 stops by `ftol` with scipy/Coleman–Li optimality 0.144307, which does not certify tight stationarity or a global optimum. The original H0 optimality is 0.210247.

## Numerical verification

The two new saved residual and Jacobian arrays replay exactly in this environment. Total and per-line objectives, nominal degrees of freedom, AICc, bounds, active-bound labels, stopping records and first-order diagnostics were independently recomputed. The summary's differences, selected paths, gate flags, line scores and continued sample exclusion all agree.

At the fixed H1 parameter vector, increasing native-pixel subdivisions from 21 to 33 and 55 changes χ² by only **+5.137 × 10⁻⁶** and **+1.103 × 10⁻⁵**, respectively. Maximum flux-model changes are 1.277 × 10⁻⁵ and 1.833 × 10⁻⁵. This quadrature effect is negligible for the reported Δχ² and adequacy decision.

Central finite differences independently checked **28 Jacobian columns**: representative log column density, velocity and log width parameters for components 0, 8 and 14; all four lines' continuum, zero-level and instrumental-width columns; and all three newly added shift columns. The largest relative L2 difference is **2.356 × 10⁻⁷**, below the declared 3 × 10⁻⁴ threshold. Derivative evaluations near active fitting bounds are mathematical checks only and do not change fitted parameters.

## Interpretation boundary

The amendment and summary correctly preserve the original failed-primary classification. The calculation reports what the allowed shifts can and cannot repair in a particular bounded gas model. It supplies neither a discovery significance nor a precision displacement, α estimate or cosmic-time point. It does not establish population-level constancy, identify the 1611 mismatch as sky contamination, or exclude other gas architectures, calibration effects, blending and continuum uncertainty. The prior numerical audit remains a historical record of the original null-only campaign; this separate review documents the subsequent authorized extension.

The machine-readable checks, source hashes, fit hashes, replay results, quadrature values and finite-difference errors are in `independent_review.json` alongside this report.
