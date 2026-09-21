# Independent validation: J233156−090802 bounded failure-case pair

Result: **PASS**, 282/282 explicit checks. No fitting was performed.

The frozen diagnostic is a bounded test of whether two regional displacements repair the previously failed conventional-profile quality screen. The original precision descriptor remains unestimated. The primary gate does not constitute a population-wide no-variation test.

## Protocol and provenance

The runner writes a non-overwritable protocol before the first additional fit. It uses one H1 fit, one H0 cross-start, and at most one H1 cross-start when H0 improves by more than 0.001. The executed branch contains two new fits; the cross-start did not improve H0 enough to trigger the optional fit. No uncertainty profile or precision datum was generated.

Frozen protocol SHA-256: `db550331c7692d802b49d957acbf5901ea02cf8ad6cc5a1b28111b2c7c8a93a5`.

Current configuration, neighboring-opacity adapter, shared engine and runner hashes match the frozen protocol. The original selected H0 JSON and NPZ still match their pre-diagnostic hashes. Additional results occupy their own directory.

## Data and model replay

The original H0 and both new fits use 26 shared gas components, 21 subpixels per native pixel, and exactly the same 1,176 globally unique native coadd pixels: 392 in each of the 1608, 2374 and 2382 regions. Independently recomputed window indices match all saved source indices; no mask change or pixel duplication occurred. Saved flux and expected-fluctuation errors match the full coadd exactly. Same-gas 1611 opacity remains inside the 1608 region without a separate likelihood.

Every saved model flux, residual, full weighted Jacobian, objective, per-line objective, nominal degrees of freedom, AICc, parameter bound and source hash was replayed. Initial parameters reproduce the documented preceding fit after the same bounded projection. H0 is exactly nested in H1 at zero regional shifts. Independent direct gradients reproduce the reported first-order diagnostics.

| Record | χ² | Nominal dof | χ²/dof | Successful stop | Explicit replay checks |
|---|---:|---:|---:|---|---:|
| n26_null_refined_os21 | 1672.841768419 | 1086 | 1.540369953 | True | 80 |
| diagnostic_H0_cross_os21 | 1672.841781033 | 1086 | 1.540369964 | True | 80 |
| diagnostic_H1_from_primary_os21 | 1667.298487118 | 1084 | 1.538098235 | True | 80 |

Selected conditional Δχ² = **5.543281301**, for two extra regional shifts. The original H0 remains the better of the two H0 endpoints. H1 still fails the same global and per-line numeric quality thresholds. This limited diagnostic does not repair the overall mismatch within the frozen architecture and ±1 km/s displacement bounds.

## Independent endpoint derivative probes

Central differences were evaluated at the actual selected H1 endpoint for the strongest weighted-Jacobian direction in each gas-parameter family, the two active Doppler-width directions, and both regional shifts. Symmetric perturbations evaluate the smooth forward model even when one side is infinitesimally outside an optimization bound; no optimizer was called.

| Parameter | Central step | Relative L2 error |
|---|---:|---:|
| logN_16 | 1e-05 | 1.4051403e-10 |
| v_21 | 0.001 | 9.684445e-09 |
| logb_16 | 1e-05 | 1.642702e-10 |
| logb_14 | 1e-05 | 1.2742384e-09 |
| shift_1608 | 0.001 | 1.8734807e-08 |
| shift_2382 | 0.001 | 1.1135475e-08 |

## Interpretation limitations

Successful ftol termination does not certify a global optimum or negligible first-order gradients. Active gas-width and zero-level boundaries remain present, and the saved gradient diagnostics remain nonzero. The 1608 nuisance moves a region containing both 1608 and 1611; it is not an isolated-transition frequency measurement. The target 2382−2374 offset remains a nonprecision diagnostic. No Gaussian discovery significance, alpha value, cosmic-age datum, or logarithmic-time fit is supported by this validation.

The finite search, chosen gas architecture, diagonal expected errors, wavelength calibration, line spread functions and blends remain limitations. The quality-screen failure itself is not evidence that every screened absorber has no physical variation.

## Exact artifact hashes

| Artifact | SHA-256 |
|---|---|
| `results/archive_complete/J233156-090802/n26_null_refined_os21.json` | `c75a19b35fdde6199296bf97dcbc95e342fde231fa6888b87d983631843d7b0b` |
| `results/archive_complete/J233156-090802/n26_null_refined_os21.npz` | `c511c9425f105020e349222b38f35595a009f8b1fb1cea187819123b59a703b2` |
| `results/archive_complete/J233156-090802/failure_pair/fits/J233156-090802/diagnostic_H0_cross_os21.json` | `6ad4bafec75d77f95f8e82a53a6ff0257295d5fd0c904f10659f1b3d8acf776c` |
| `results/archive_complete/J233156-090802/failure_pair/fits/J233156-090802/diagnostic_H0_cross_os21.npz` | `72a999498563a091a71cc3b81cd2227057d12123572f7888e53bcf6a8ad91289` |
| `results/archive_complete/J233156-090802/failure_pair/fits/J233156-090802/diagnostic_H1_from_primary_os21.json` | `f0ba0c57434c769d993dfa94cbb9196904e577007f6634c78cb1b18d4f7760ac` |
| `results/archive_complete/J233156-090802/failure_pair/fits/J233156-090802/diagnostic_H1_from_primary_os21.npz` | `1fae329e9cb62ddb0290edb61644d3f4a7fd43278ce26d3ace68f9f994ddd899` |

## Failed checks

None.
