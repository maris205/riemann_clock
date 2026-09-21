# J064326−504112 independent numerical fit review

**Result: 304 / 304 grouped mechanical checks passed across all 20 saved fits.** This validates reproducibility of the saved calculation, including its deliberately negative adequacy decision. It does not certify a stationary/global optimum, identify the origin of the weak-line residuals, or establish a change in any physical constant.

The review used `archive_complete_model.Model` without modifying that adapter, its underlying Voigt implementation, input sources, configurations, or saved fit products. No optimizer was run. The checked scientific quantities and whole-array comparisons are enumerated individually in `results/archive_complete/J064326-504112/independent_fit_review.json`; individual pixels are not counted as separate tests.

## Scope and source consistency

The directory contains exactly 18 initial ordinary-model fits: component counts 4, 6, 8, 10, 12 and 16, each with starting seeds 0, 1 and 2. These are followed by the selected 16-component ordinary fit at 21 subdivisions per native pixel and one separate ordinary-only control omitting Fe II 1611. All 20 records have `free_shifts=false`, no shift parameter labels, empty shift dictionaries, and no fixed-shift profile. No free-shift fit NPZ or selected H0/H1 pair exists.

For every saved fit, the current coadd, metadata, atomic table, base-model source and adapter hashes match the recorded hashes. The audit independently reconstructs the isotope-weighted wavelength references and the fixed velocity selection from the full coadd. Every retained source index, valid mask and expected-fluctuation error agrees. Saved wavelengths, velocities, fluxes, both error arrays, masks, source indices and predicted profiles also agree with the reconstructed model.

The primary retains exactly 436 native pixels: 109 for each of Fe II 1608, 1611, 2374 and 2382. All selected source-valid expected errors remain positive and finite. No fitted-residual masking occurs.

At every saved parameter vector, direct replay reproduces the saved residual and analytic Jacobian arrays exactly in this environment. The total objective, each line's objective, nominal parameter count and degrees of freedom, AICc, parameter bounds, active-bound labels and recorded first-order diagnostics were recomputed. The saved optimizer stops are successful status-2 (`ftol`) stops; this is checked as a record of termination, not treated as proof of tight stationarity.

## Selection and adequacy decision

The chosen start at every component count is the minimum-χ² successful ordinary-model start. Selection across counts uses ordinary-model AICc only. The lowest tested AICc occurs at 16 components; the refined primary is `n16_null_refined_os21`.

| Primary quantity | Recomputed value |
|---|---:|
| Native pixels / parameters | 436 / 64 |
| χ² / nominal degrees of freedom | 532.225273754 / 372 |
| χ² per nominal degree of freedom | 1.430713101 |
| Fe II 1608 χ² per pixel | 0.865216871 |
| Fe II 1611 χ² per pixel | **2.058481575** |
| Fe II 2374 χ² per pixel | 1.074081611 |
| Fe II 2382 χ² per pixel | 0.885020620 |

The declared conditional gate requires successful stopping, total χ²/ν ≤ 1.5 and every line's χ²/pixel ≤ 1.8. The primary passes the aggregate threshold and fails the line threshold solely for 1611. This means that **1611 alone triggers this numerical gate**; it does not establish that 1611 has a unique physical failure mechanism. The saved failed-gate decision and withholding of relative-shift fits are correct.

Two caveats remain important. The selected primary has Coleman–Li/scipy optimality 0.210247414, normalized gradient mapping 0.351857254 and column-normalized gradient maximum 0.613461489. It is not a gradient-tolerance solution. `logb_14` is active near the b = 0.5 km/s lower bound, while `zero_1611` is active at +0.02. Moreover, 16 is the largest searched component count, with AICc still declining over the grid. The result describes a bounded exploratory architecture search, not an exhaustive determination of the absorber's gas structure. Different starts at counts 8–16 also leave different local minima.

## Numerical quadrature and derivatives

The selected 21-subdivision parameter vector was evaluated unchanged at 33 and 55 subdivisions. These are quadrature checks, not additional fits.

| Subdivisions | Fixed-parameter χ² | Change from 21 subdivisions | Maximum flux-model change |
|---:|---:|---:|---:|
| 21 | 532.225273754 | 0 | 0 |
| 33 | 532.225279073 | +0.000005319 | 0.000012859 |
| 55 | 532.225284990 | +0.000011236 | 0.000018460 |

This numerical integration difference is immaterial to the reported gate failure. It does not remove physical-model or optimization uncertainty.

Independent central finite differences checked 25 representative residual-Jacobian columns: log column density, velocity and log width for components 0, 8 and 14, together with continuum constant, continuum slope, zero level and log instrumental width for every line. The largest relative L2 derivative difference was **4.185 × 10⁻⁶**, below the stated 3 × 10⁻⁴ tolerance. For active-bound parameters, a small excursion across the fitting bound was used solely to check the mathematical derivative; no fitted parameter was changed.

## Ordinary-only omission control

`omit1611_n16_null_only_os21` retains exactly the same 327 pixels, fluxes, uncertainties, masks and source indices for 1608, 2374 and 2382, with the same 16 gas components and common velocity window. Its initial gas parameters match the primary after the documented clipping infinitesimally inside the bounds. It has no additional relative-shift parameters and does not reselect the component count.

The control gives χ² = **300.248362083** with 267 nominal degrees of freedom and 60 fitted parameters. Its three χ²/pixel values are 0.895065667, 1.056396796 and 0.803109666, so the same conditional gate passes. There are no recorded active bounds. Its optimality is 0.059592921 with `ftol` termination, again not a proof of a global optimum.

This post-primary diagnostic supports the narrow statement that the remaining three lines admit a substantially better ordinary fit under the same selected component count. It also removes the useful unsaturated 1611 constraint. It **does not replace the failed four-line primary**, identify the 1611 residuals as sky contamination, or authorize a new-physics claim. The 1611 sky-proxy warning was recorded before fitting and remains a reason to investigate its data quality, not a confirmed explanation.

The errors throughout remain diagonal expected-fluctuation weights. Pixel covariance, exposure-specific sky and atmospheric response, calibration distortions, continuum choices and gas architecture uncertainty remain outside this mechanical audit. No α or cosmic-time-law measurement has been made.
