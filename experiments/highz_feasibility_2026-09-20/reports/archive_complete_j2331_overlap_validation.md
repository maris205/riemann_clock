# Independent pre-fit validation of J233156−090802 neighboring opacity

Date: 2026-09-21. Scope: the approved [−480, +500] km/s selection and target-specific `archive_complete_j2331_model.py` adapter. No optimization was run, no selected fit was inspected, and no engine or selection source was edited by this review.

**Result: PASS for the tested adapter, native-pixel selection, and gas/region-shift derivatives.** The known same-gas 1611 overlap is now included in the 1608 spectral region. This resolves the omission identified in `archive_complete_j2331_selection_review.md`; it does not establish that all other blends or instrumental effects are absent.

## Frozen inputs checked

| File | SHA-256 at validation |
|---|---|
| `results/archive_complete/J233156-090802/selection.json` | `5202efc55d480498f70230755819284699e9bf8685caf48ae0c2cac55742aaeb` |
| `code/archive_complete_j2331_model.py` | `9167a0fffb5bbb1e05b252241e23a05f56d72c11a133182b67a9904bbc5c2a20` |
| `code/archive_complete_model.py` | `6bd3871b82bc61c4fea1d8e8c91bc22db78c055c06770d0cc242ef39f5182416` |
| `code/archive_profile_pilot.py` | `f6b0f8ef92d5bc5e00a02dee2d4ce6e386836af78712148823d9deba5992ce20` |
| `data/raw/espresso_null/MM_VPFIT_2013-11-10.dat` | `0bbf88665e1e6d108fc9390d7854ccb1fdb97e993cb9ca80183b65b27b52983a` |

The full processed-coadd hash still equals the selection record: `80276c58831dc457e9c7f02521ef4883311fa9c2057e5a4306454c607c13b4d0`. The current selection-script hash also equals the recorded value. The earlier independent raw tar/FITS/processed-array checks remain documented in the selection review.

## Adapter and native-pixel checks

An independent script performed **71 explicit Boolean checks, all passing**. These covered recorded hashes, class/evaluator identity, native-pixel selection, all preserved region attributes, atomic rows and likelihood uniqueness.

The adapter appends precisely the four Fe II 1611 isotope rows independently parsed from the atomic text file to the four existing 1608 rows. The atom arrays for 2374 and 2382 are unchanged. The inherited evaluator uses the same gas N, v and b for these rows and accumulates their optical depths before computing transmission, Gaussian convolution and finite-pixel integration. It does not add independent 1611 gas parameters or a fourth likelihood region.

Comparison with the unmodified campaign model confirmed exact equality of each region's reference wavelength, velocity samples, dispersion, oversampled integration grid, convolution padding, source-valid mask, flux, expected-fluctuation errors, statistical errors, wavelengths, global source indices, continuum coordinate and nominal/arc resolution metadata. The three likelihood regions remain 1608, 2374 and 2382.

Using an independent logarithmic-velocity calculation from the full coadd, the inclusive expanded interval selects:

| Region | Global source indices, zero based and inclusive | Retained native pixels |
|---|---:|---:|
| 1608, with neighboring 1611 opacity | 51590–51981 | 392 |
| 2374 | 98298–98689 | 392 |
| 2382 | 98717–99108 | 392 |

All **1,176 global source indices are unique**. H0 and H1 use exactly the same indices and validity masks. The data values and expected errors equal the corresponding full-coadd values exactly. No likelihood pixels are duplicated by the inclusion of neighboring opacity.

For the six-component derivative probe, H0 has 30 parameters and H1 has 32, as in the unmodified architecture. The existing 1608-region nuisance shift moves both 1608 and 1611 opacity within that observed region. This is explicitly a spectral-region displacement; it is not a separately measured 1608 transition frequency or alpha parameter. The 2382-versus-2374 diagnostic retains its original parameter definition.

## Independent finite-difference checks

The probe used 21 subpixels per native pixel and six gas components, at a deliberately specified parameter point rather than a fitted solution:

- log10 N: `[14.7, 14.2, 14.8, 14.4, 13.5, 13.3]`.
- Gas velocities in km/s: `[−400, −320, −140, −80, 72, 360]`.
- Doppler b in km/s: `[10, 12, 15, 8, 6, 8]`.
- Region shifts in km/s: 1608 = `+0.24`, 2382 = `−0.18`.
- Other parameters: the model's seed-zero initial nuisance values.

All 18 gas columns and both existing region-shift columns were compared with symmetric central differences of the weighted residual. Steps were 0.001 km/s for gas velocity and region shift, and 0.00001 for log10 N and log b. The relative discrepancy was computed as `norm(J_FD − J_analytic) / norm(J_analytic)` across all likelihood pixels.

| Parameter family | Number of checked columns | Largest relative L2 discrepancy |
|---|---:|---:|
| log10 N | 6 | 1.56 × 10⁻¹⁰ |
| Gas velocity | 6 | 8.74 × 10⁻⁹ |
| log b | 6 | 3.16 × 10⁻¹⁰ |
| 1608-region shift | 1 | 5.95 × 10⁻⁹ |
| 2382-region shift | 1 | 7.88 × 10⁻⁹ |

To exercise the neighboring contribution specifically, the same parameters were evaluated with the isolated-1608 base model. For each checked direction, the analytic difference between neighbor-model and base-model Jacobian columns was compared with a central difference of their residual difference. Every derivative-difference error had absolute L2 norm below **5.1 × 10⁻⁹**. The largest relative discrepancy, 3.46 × 10⁻⁵, occurs for a nearly vanishing neighboring contribution with norm only 1.39 × 10⁻⁴; the full derivative remains accurate to approximately 1.5 × 10⁻¹⁰. The 2382 shift has exactly zero neighboring contribution in the 1608 residual block, as required.

These are finite-point numerical derivative checks, not evidence of optimization convergence, model identifiability or a global optimum.

## Runner behavior checked

Importing the adapter alone does not replace the campaign engine's Model class. The runner applies its replacement within the runner process. A separate no-fit process-pool probe reproduced this sequence under the current environment's `fork` start method: the child constructed `NeighborModel`, contained eight isotope rows in the 1608 region, and retained the same 1,176 likelihood pixels. The probe restored the class in its own process afterward; no source file was modified.

This verifies the current Linux/fork execution path. A port to a process-start method that reloads modules independently, such as `spawn`, would need an explicit adapter initializer or equivalent child configuration; the inheritance of an in-memory class replacement should not be assumed there.

## Scope of the pass

The reviewed implementation supplies the known neighboring opacity without changing observations, weights, primary reference frames, or parameter counts. The fixed expanded selection addresses the previously omitted detached red candidates. Remaining inference still depends on gas architecture, saturation, continuum, instrument response, calibration and other blends. This pre-fit pass authorizes no scientific claim about measured shifts, alpha, or a logarithmic cosmic-time law.
