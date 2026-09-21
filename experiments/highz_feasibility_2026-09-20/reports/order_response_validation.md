# Centered order-response model: independent numerical validation

## Material Passport

- Mode: reproducibility and statistical interpretation validation.
- Material: public native ESPRESSO HE 0515−4414 spectra, fixed shared Fe II gas template, and newly generated conditional order-response fits.
- Review provenance: separate agent in the same workspace and conversation context; an additional agent independently derived the kernel moments. This is not an external or blinded audit.
- Verification status: **VERIFIED within the stated computational scope**.
- Checks: 7278 / 7278 passed; 0 failed.
- Review artifacts: `code/order_response_validate.py` and `results/order_response/validation.json`.

The audit reconstructs the response by direct sampled-kernel convolution and independently solves each row's linear nuisance parameters with SVD least squares. Production uses shifted Gaussian-filter splines and QR profiling. Saved objectives, residuals, source bindings, matched-data nesting, parameter counts, covariance calculations, and chronological response transfer were checked. This validates the specified conditional calculation, not the physical adequacy of its noise, gas, or instrument assumptions.

## Kernel and identifiability

Writing s = F / 2.354820045 and d = 1.2 cbrt(a), the 80:20 mixture has area 1, mean 0, variance s², third moment 0.165888 a, and fourth moment 3s⁴ + 0.01327104 |a|^(4/3). Direct continuous quadrature verifies these identities. The [1.6, 2.6] km/s equivalent-width bounds keep the component variance strictly positive; the limiting requirement is F > 1.1303136216 km/s.

F denotes the Gaussian-equivalent width determined by the second moment, not the actual asymmetric kernel FWHM. The parameter a controls the third moment and is not standardized skewness. The complete mixture is C¹ but generally not C² at a = 0. Its first derivative there is −0.027648 G‴. Checked finite differences are appropriate for this numerical diagnostic; a naive derivative through cbrt(0) is not.

Zero kernel centroid removes an exact kernel translation, but does not make fitted shift and asymmetry independent. Even isolated Gaussian absorption has a shift/skew tangent correlation of magnitude √(3/5) ≈ 0.775 under uniform continuous noise. Finite windows, blends, saturation, and nuisance profiling can alter that coupling. An asymmetric LSF can move a line core while leaving its full unmasked centroid unchanged.

## Reproduction results

- Main nonlinear fits audited: 40.
- Fixed-offset nonlinear profile fits audited: 52.
- Maximum direct-kernel standardized-residual discrepancy: 2.18868e-08.
- Maximum direct-kernel objective discrepancy: 7.01744e-07.
- Full local Jacobian rank and a separate covariance inverse were checked for every main fit; both selected full A1 fits also received all-column centered finite-difference checks.
- Refinement used half the velocity-grid spacing, 21-point pixel quadrature, and a separately differenced optimizer Jacobian for actual nonlinear re-fits.

| Transition | Refined − saved offset (m/s) | Refined − saved χ² | Refit success |
|---|---:|---:|---|
| 2600 | -0.00168117 | -8.73639e-07 | True |
| 2382 | -0.00782136 | -7.73311e-07 | True |

Multiple starts are retained; selection uses the minimum objective among successful starts. The largest asymmetry-family main-fit spread was 5.34419e-06 in χ². Optimizer termination, numerical reproducibility, and agreement among these starts do not prove a global minimum.

## Conditional interpretation

For Fe II 2600, the selected Gaussian response gives an order difference −62.75 ± 16.66 m/s; adding the centered asymmetric response gives −58.64 ± 18.65 m/s. These are conditional local standard errors under diagonal native ERRDATA and a fixed gas template. The zero-offset versus free-offset improvement decreases from 14.21 to 10.86 in χ². Thus this particular response family does not erase the order difference. It also does not establish that the residual difference is beyond conventional instrument, extraction, wavelength calibration, or gas-model effects.

The full 2600 A0 fit and several training A0 starts touch an asymmetry bound. This restricts interpretation of its conditional profile; shape bounds are assumptions rather than measured calibration limits. The Fe II 2382 control remains much closer to zero and has a weak order-offset preference. Full minimum reduced objectives are about 1.34 (2600) and 1.29 (2382), so model/noise adequacy is not established by these fits.

The local A1 shape estimates are weak: for 2600, a52 = −0.353 ± 0.516 and a53 = −0.925 ± 0.775; for 2382, a41 = −0.341 ± 0.990 and a42 = −0.044 ± 1.637. These conditional local errors ignore the hard [−1, 1] bounds and are not calibrated confidence intervals. In particular, an inactive bound flag does not imply that the nearby boundary is inferentially irrelevant. The shapes are not empirical LSF measurements.

Training on 2018 exposures and freezing response parameters on 2019–2020 is verified exactly, with disjoint exposure indices. For 2600, A1 yields held-out χ² larger by 5.15 than G1. This is a conditional response-transfer check, not independent prediction: the fixed gas template was estimated from the full coadd, including these held-out photons. It also supplies no second cosmological-age sample.

An astrophysical frequency displacement of one transition is common to its simultaneous measurements in both orders. A difference between orders in the same exposure is therefore a useful conventional-effect diagnostic, but this sensitivity experiment does not isolate which conventional mechanism produces it. No cosmic-time-law or varying-constant inference is calibrated here.

## Statistical fallacy scan — 11 / 11 considered

| Item | Scope-specific assessment |
|---|---|
| Simpson's paradox | Shared response is an aggregate assumption. Prior exposure/chronological contrasts and the held-out loss are retained; pooled results do not establish identical exposure-level response. |
| Ecological fallacy | No inference from a pooled order diagnostic to individual atomic evolution or multiple cosmological ages. |
| Berkson/selection bias | This target and transition were chosen after the earlier order discrepancy; the diagnostic is conditional on that selection. |
| Collider bias | No causal adjustment identification is claimed; nuisance profiling cannot establish a causal LSF mechanism. |
| Base-rate neglect | Not a diagnostic classifier; no posterior discovery probability inferred from an objective difference. |
| Regression to the mean | Reusing the selected discrepancy motivates a control transition and retained chronological transfer results; shared-template reuse remains explicit. |
| Survivorship bias | All planned main starts and grid points are retained, including boundary fits; the minimum-successful rule and total counts are audited. |
| Look-elsewhere effect | Exploratory model family and previous searches preclude treating local objective improvement as a calibrated discovery probability. |
| Forking paths | Protocol frozen before these model comparisons, but after seeing earlier data. Assumed mixture weights/amplitude and all four models are disclosed. |
| Correlation versus causation | Better fit or same-order direction does not identify physical or instrumental cause. |
| Reverse causality | No causal or temporal ordering is inferred; observing years are not different source cosmic ages. |

## Re-run

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_clock/experiments/highz_feasibility_2026-09-20/code/order_response_validate.py
```

The JSON records exact source/result hashes, checks, numerical refinement metrics, active bounds, and multistart spreads. Prior raw FITS validation is linked by artifact hashes rather than represented as a new raw download or fresh independent astrophysical dataset.
