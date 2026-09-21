# J064326−504112: completed bounded conventional-model campaign

**Outcome: the four-line primary fit does not pass its predeclared per-line adequacy gate. No relative-shift alternative, alpha estimate or cosmic-time point is reported.**

The primary line set (1608, 1611, 2374, 2382), common [−145,+100] km/s window and initial gas-center maps were frozen before fitting shifts. All 436 source-valid native pixels were retained. Errors are the original FITS expected-fluctuation row with diagonal weighting. Weak 1611 was included from the outset to constrain saturation, with its pre-fit 5897 Å strong-sky proxy warning preserved.

| Gas components | Best ordinary χ² | Nominal dof | AICc | Successful stop |
|---:|---:|---:|---:|:---:|
| 4 | 4294.841839 | 408 | 4354.832011 | True |
| 6 | 1526.090181 | 402 | 1600.025343 | True |
| 8 | 879.908839 | 396 | 968.212637 | True |
| 10 | 623.023617 | 390 | 726.139298 | True |
| 12 | 581.701547 | 384 | 700.093192 | True |
| 16 | 532.227016 | 372 | 682.652892 | True |

Three frozen starts were run per count (4, 6, 8, 10, 12, 16), with at most 350 objective evaluations each. All 18 initial fits stopped successfully; minimum ordinary-model AICc selected 16 components, the maximum searched count. Refinement at 21 subpixel integration points gives χ² = 532.225273754 for 372 nominal degrees of freedom. This does not establish an optimum over all gas architectures.

| Transition | Pixels | Ordinary χ² | χ² / pixel |
|---|---:|---:|---:|
| 1608 | 109 | 94.308639 | 0.865217 |
| 1611 | 109 | 224.374492 | 2.058482 |
| 2374 | 109 | 117.074896 | 1.074082 |
| 2382 | 109 | 96.467248 | 0.885021 |

The predeclared diagnostic gate requires successful stopping, total χ² / nominal dof ≤ 1.5, and every transition χ² / pixel ≤ 1.8. The weak 1611 line fails the last condition. The high-column core is saturated in the other transitions, so an apparently acceptable fit to those lines does not resolve the residual mismatch in the weak constraint. Possible sky residuals, unrelated absorption, continuum treatment, gas architecture or line response require further examination. This is not evidence that a varying constant is required.

The initial selection audit independently reproduces all 436 native pixels from raw FITS. The selected blue boundary is a defensible local separation point, but nearby features are not identified and the main complex is not claimed to be an isolated complete absorber. Edge residual diagnostics are retained without clipping or changing the frozen window.

A null-only explanatory omission of weak 1611, if present, must remain a separate control. It cannot overwrite the failed four-line primary fit or authorize a physical variation interpretation.

Files: `config.json`, `campaign_protocol.json`, `ordinary_model_summary.json`, all individual fit JSON/NPZ files, `ordinary_profile_fit.pdf`, `ordinary_count_comparison.pdf`, and `edge_residual_diagnostics.json` under `results/archive_complete/J064326-504112/`.

The largest weak 1611 residual is at 5893.306834 Å (velocity -105.553 km/s), 5.017 in the saved standardized residual. It is outside the conservative sky proxy. The proxy flag alone therefore does not establish an explanation for the mismatch. `weak1611_wavelength_diagnostic.pdf` shows the wavelength mapping; all pixels remain in the primary fit.

The selected primary has active `logb_14` (lower 0.5 km/s) and `zero_1611` (upper +0.02) bounds. Its recorded SciPy first-order optimality is 0.210247, so successful ftol termination is not claimed as strong stationarity or a certified global minimum.

The additional null-only weak 1611 omission retains the fixed 16 gas architecture and the other 327 native pixels. It gives χ²=300.248362083 for 267 nominal degrees of freedom, successful stopping=True, and the same conditional adequacy gate=True. This localizes sensitivity to including the weak constraint, while also losing saturation information; it does not identify contamination and does not replace the failed primary. No displacement alternative was fitted.

Independent validation: 304/304 grouped mechanical checks PASS across all 20 saved fits, including exact objective/residual/Jacobian replay, raw source-pixel/error identity, and 25 representative finite-difference columns. At fixed selected parameters, increasing integration from 21 to 55 subpixels changes χ² by only 1.12e-05. These checks validate calculation and reproducibility, not model completeness or physical variation. See `archive_complete_j0643_fit_review.md`.
