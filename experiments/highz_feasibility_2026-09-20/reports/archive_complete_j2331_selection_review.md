# Independent selection review: J233156−090802

Reviewed on 2026-09-21. Scope: read-only inspection of the frozen selection, the full coadded source, the ±500 km/s inspection image, and the earlier apparent-optical-depth (AOD) image. No profiles were fitted, no shift results were inspected, and no data or selection files were changed by this review.

**Disposition after the pre-fit addendum below: source and pixel checks pass, but a physical model that treats the expanded 1608 region as an isolated transition is inadequate. The same gas produces overlapping Fe II 1611 opacity. Even the original red margin can be affected. Model the neighboring transition or explicitly justify and freeze a conservative per-window exclusion before fitting.**

The initial dominant-complex assessment below is preserved as the audit history. Its interpretation of 1608 red depressions and continuum margins is qualified by the subsequently discovered neighboring-transition overlap. Neither version verifies selection of every Fe II component belonging to the absorber.

Subsequent pre-fit resolution: the approved window was expanded to [−480, +500] km/s and a same-gas 1611 adapter was added. Its independent source, pixel and derivative validation passed; see `archive_complete_j2331_overlap_validation.md`. The hashes and 756-pixel counts below refer to the earlier archived selection, whereas the validated expanded selection contains 1,176 unique pixels.

## Source and pixel checks

The source named by `results/archive_complete/J233156-090802/selection.json` is the full processed coadd, `data/processed/archive_expansion/J233156-090802.npz`, rather than one of the earlier ±150 km/s extracts. Its SHA-256 is:

`80276c58831dc457e9c7f02521ef4883311fa9c2057e5a4306454c607c13b4d0`

This matches the selection record. Independent extraction of the FITS bytes from the downloaded tar archive reproduced the metadata hashes:

- Tar archive: `07e4ae515bb8356cc89c843e2b0bc17a616d4ff8cbd3608424ef191d5da306b6`.
- FITS member: `034bf70596baa5ff0079ce8c878875fc65bb36bc01175a67ddc1c57470f2cad3`.

The processed flux, statistical-error, expected-fluctuation, and status arrays are exactly equal to the corresponding original FITS arrays. Independently reconstructing the wavelength array from CRVAL1, CRPIX1 and CD1_1 gives zero difference from the processed array.

The velocity coordinate was recomputed as `299792.458 * ln(wavelength / [weighted_atomic_reference * 3.143])`, using the oscillator-strength-weighted isotope reference for each transition. The common inclusive interval [−480, +150] km/s yields:

| Transition | Native coadd indices, inclusive and zero based | Pixels inside interval | Valid pixels with finite flux and positive finite expected fluctuation |
|---|---:|---:|---:|
| Fe II 1608 | 51590–51841 | 252 | 252 |
| Fe II 2374 | 98298–98549 | 252 | 252 |
| Fe II 2382 | 98717–98968 | 252 | 252 |

Thus the specified selection retains all **756 native pixels** across the three lines. There are no internal invalid pixels or residual-based exclusions in this selection definition. This review checks the selection itself; it does not certify a subsequently saved fit mask.

At review completion, the selection JSON SHA-256 was `91719cd1cd2fdc4f4742aeac89634b655e10f507eb69f9a0410ffe7da119c1c9`. Its recorded selection-script SHA-256, `215c307f0249d288da7d19a01b85d6604114d7aa2d52a5ca017f10322d295247`, matched the then-current script. Later cosmetic updates to the plot title may change these hashes and should be recorded normally.

## Scientific assessment of the window

The ±500 km/s image shows a broad, recognizably shared Fe II complex extending approximately from −430 to +85 km/s, with the strong 2382 transition substantially saturated over much of that range. The [−480, +150] km/s interval covers that dominant complex and includes useful outer pixels on each side. The weaker 1608 and 2374 transitions retain considerable shape information where 2382 is saturated. Joint conventional-gas modelling is therefore a reasonable exploratory next step.

The earlier AOD image only covers approximately [−150, +150] km/s. It supports shared absorption within that limited range, but cannot establish the extent of the blue part of the full complex. Its suppressed 2382 AOD is expected under saturation and the stated low-flux floor; it is not evidence for a frequency anomaly.

The 1608 continuum is visibly below one. The selection acknowledges this and permits continuum flexibility. The outer margins should not be described as a precise independent continuum determination; line-dependent continuum error and curvature remain potential limitations. The large saturated extent also makes component decomposition non-unique. A finite count grid up to 26 components is a documented search scope, not proof of gas-model completeness.

## Detached red features: important qualification

The troughs near +180, +360 and +450 km/s are strongest in 2382, but the data do **not** establish that they are absent from Fe II 2374 or 1608. To check the visual impression without any nonlinear fitting, this review calculated inverse-variance weighted raw flux means using the expected-fluctuation array:

| Velocity interval (km/s) | 1608 mean flux | 2374 mean flux | 2382 mean flux |
|---|---:|---:|---:|
| +160 to +205 | 0.85885 ± 0.01241 | 0.93527 ± 0.01015 | 0.61286 ± 0.00895 |
| +345 to +385 | 0.79799 ± 0.01342 | 0.92059 ± 0.01270 | 0.64289 ± 0.00907 |
| +440 to +470 | 0.88273 ± 0.01656 | 0.94523 ± 0.01306 | 0.73217 ± 0.01091 |
| +210 to +320, comparison interval | 0.88672 ± 0.00824 | 0.97823 ± 0.00680 | 0.96860 ± 0.00665 |

These errors assume independent pixels and omit continuum uncertainty; they are descriptive diagnostics, not detection significances. The 2374 depressions in the three detached intervals and some corresponding 1608 structure preclude calling the troughs “2382-only.” The summed atomic fλ ratio is approximately 10.26 for 2382/2374, so a large 2382 absorption depth can coexist with a relatively small 2374 depression. Such behaviour is compatible with weak detached Fe absorption and does not require an external blend identification.

The reviewed selection was amended before fitting to call the interval the **dominant visible shared complex** and explicitly leave detached red features as weak-Fe-or-external candidates of unidentified origin. That wording is supported. Their exclusion is acceptable as a fixed, disclosed scope restriction for exploration of the dominant complex. It is not evidence that the excluded features are unrelated. A wider-window sensitivity analysis, or a justified identification of those features, would be necessary before any claim that the fitted model accounts for the complete absorber.

The inspection PNG initially retained “complete-system selection” in its title even after the JSON wording was amended. The fitting agent was asked to change this to “dominant-complex selection.” This is a presentation correction; the numerical window remains unchanged.

## Interpretation boundary

This is a data-inspected, exploratory selection, not a blind confirmatory sample. An adequate fit would establish conditional compatibility of the selected dominant-complex pixels with a conventional shared-gas model. An inadequate fit would first require investigation of gas complexity, continuum, line spread functions and blends. Neither outcome alone provides a measurement of alpha, a cosmic-age dependence, or a test of the proposed logarithmic time law.

## Pre-fit addendum: Fe II 1608 and 1611 overlap

The fitting agent identified a neighboring-transition issue before performing the proposed expanded-window fits and requested an independent calculation. Re-reading all isotope entries from the recorded atomic file gives oscillator-strength-weighted vacuum reference wavelengths:

- Fe II 1608: 1608.4508522723688 Å, with summed f = 0.05769999.
- Fe II 1611: 1611.200368969589 Å, with summed f = 0.00138.

For the logarithmic velocity coordinate used in this analysis, the separation calculated from these weighted references is

\[
\Delta v = c\ln(\lambda_{1611}/\lambda_{1608})
         = 512.0334549977579\ {\rm km\,s^{-1}}.
\]

Thus a gas component at velocity u in the 1611 frame contributes at u + 512.033455 km/s in the 1608 frame. The common redshift cancels from this mapping. Individual isotope offsets remain part of the physical opacity calculation; the value above is the exact mapping between the two specified weighted reference frames.

| Region or feature | Corresponding interval in the other frame (km/s) |
|---|---:|
| Dominant gas in 1611, u = −430 to +85 | 1608 frame: +82.033455 to +597.033455 |
| Proposed 1608 window [−480, +500] | 1611 frame: −992.033455 to −12.033455 |
| Original 1608 window [−480, +150] | 1611 frame: −992.033455 to −362.033455 |
| 1608 interval +115 to +150 | 1611 frame: −397.033455 to −362.033455 |
| 1608 interval +160 to +205 | 1611 frame: −352.033455 to −307.033455 |
| 1608 comparison interval +210 to +320 | 1611 frame: −302.033455 to −192.033455 |
| 1608 interval +345 to +385 | 1611 frame: −167.033455 to −127.033455 |
| 1608 interval +440 to +470 | 1611 frame: −72.033455 to −42.033455 |

Consequently, the expanded 1608 window intersects much of the dominant absorber's 1611 profile. Some of the earlier 1608 red depressions can be same-gas 1611 absorption; they are not clean evidence for additional components in the 1608 transition at those plotted velocities. Likewise, the mean 1608 flux in +210 to +320 km/s is not a clean continuum measurement because that interval coincides with 1611 velocities −302 to −192 km/s. The independent depressions in 2374 remain a reason not to label the detached 2382 troughs “2382-only,” but the 1608 columns in the earlier table must now be interpreted with this overlap caveat.

There is a reciprocal issue in the inspection panel labelled 1611: features plotted around −490 and −450 km/s correspond to 1608 velocities +22.033 and +62.033 km/s. Hence the blue absorption in that panel cannot be established as an unrelated external blend; overlapping 1608 from the same absorber is a conventional explanation that must be checked first. The excluded-line reason should avoid asserting an unrelated identification.

The summed fλ ratios are 0.0239576991 for 1611/1608 and 0.0299170958 for 1611/2374. These are small, but they do not make 1611 negligible in a strongly saturated, high-column-density complex. Its opacity is physically specified by the same gas parameters rather than an arbitrary continuum correction.

**Preferred treatment:** retain the predeclared expanded observed-pixel window and add all 1611 isotope components to the 1608 region's opacity using the same Fe II column densities, gas velocities and Doppler parameters. Sum the 1608 and 1611 opacities before exponentiation, instrumental convolution and pixel integration. No separate gas-column freedom is needed for this same-ion neighboring transition. Any later diagnostic shift parameter must explicitly state whether it shifts an observed wavelength window or one transition; those are distinct when a window contains two transitions. If a separate 1611 pixel window is also added, overlapping native pixels must not be counted twice.

**Alternative restricted treatment:** freeze a shorter 1608 window and retain the wider 2374 and 2382 windows, with the restriction disclosed as line-dependent coverage. A cutoff at +80 km/s is not a proof of zero contamination. The permitted gas center −465 km/s maps to +47.033455 km/s in the 1608 frame; finite Doppler and instrumental wings extend farther blueward. Even for visible centers beginning around −430 km/s, line wings extend blueward of +82.033455. A safe cutoff therefore requires a stated wing/error tolerance and justified gas bounds, not just the nominal center separation. Using the known neighboring opacity is the cleaner default for the proposed expansion.

These findings were sent directly to both the fitting agent and `/root/clock_archive_expansion` before fitting. This addendum performs no optimization and establishes no measured relative displacement.
