# Independent validation of the pixel-noise pilot

**PASS: 63 checks.** No material numerical error was found. This review verifies a synthetic-line benchmark, not an astronomical measurement.

The five analytic Gaussian derivatives agree with independent five-point finite differences in all twelve data windows and two additional profile configurations. All twelve window hashes match the benchmark manifest. An independent Cholesky whitening and nuisance-space projection reproduces the full inverse-Fisher result, including depth, width, continuum offset and continuum slope. The adjacent-pixel covariance is positive definite.

| Target / mock transition | Independent diagonal σ, m/s | With illustrative adjacent covariance, m/s | MC run |
|---|---:|---:|---|
| J051707-441055 / FeII_1608.4509 | 72.373 | 89.685 | Yes |
| J051707-441055 / FeII_2344.2128 | 40.251 | 49.885 | Yes |
| J051707-441055 / FeII_2374.4601 | 40.344 | 49.995 | Yes |
| J051707-441055 / FeII_2382.7642 | 38.739 | 48.006 | Yes |
| J051707-441055 / FeII_2586.6493 | 37.264 | 46.179 | Yes |
| J051707-441055 / FeII_2600.1725 | 39.136 | 48.503 | Yes |
| J034943-381030 / FeII_1608.4509 | 107.497 | 133.211 | Yes |
| J034943-381030 / FeII_2344.2128 | 217.491 | 269.518 | Yes |
| J034943-381030 / FeII_2374.4601 | 371.387 | 460.278 | Yes |
| J034943-381030 / FeII_2382.7642 | 431.556 | 534.847 | Yes |
| J034943-381030 / FeII_2586.6493 | 7055.714 | 8743.528 | No; low SNR |
| J034943-381030 / FeII_2600.1725 | 4480.401 | 5552.762 | No; low SNR |

The illustrative correlation increases the centroid uncertainties by about 24%. All ten eligible Monte Carlo summaries are compatible with their stated local Gaussian benchmarks at the precision supported by 300 trials. No failed or boundary-constrained fit was recorded. The archived MC summaries were reviewed; the nonlinear MC was not rerun in this independent short audit.

## Interpretation requirements

- This calculation fits only synthetic profiles, never the archived observed flux. Its outputs are local statistical precision benchmarks for the specified mock line, not measured astronomical line errors or constraints.
- The assumed R=45,000, intrinsic b=3 km/s and post-instrument depth 0.25 are fixed choices. Archived coadds may have varying instrumental profiles, and actual lines can be weaker, saturated, blended or multicomponent.
- The mock noise is the median normalized error over the full archived window and is then held constant. This approximates a usable noise scale; it is not a refit of actual pixel-dependent photon statistics, continuum uncertainty, or pipeline error calibration.
- Correlated noise is an illustrative positive-definite tridiagonal covariance with adjacent-index rho=0.3. It is not an estimate of the archive covariance, and the nonlinear MC was run only for independent noise.
- All four specified profile nuisances are marginalized in the Fisher calculation and fitted in the nonlinear Monte Carlo. Additional real calibration and atomic uncertainties are excluded.
- Each MC has 300 trials: coverage estimates fluctuate by about 1.26 percentage points under a true 95% coverage model. Individual observed coverages of 0.933–0.960 do not indicate a resolved failure or establish precise tail probabilities.
- The two z≈3.025 windows with SNR 1.31 and 2.06 are correctly excluded from MC. Their 4.5–7.1 km/s local-Fisher values must not be treated as attainable centroid errors for detections.
- The high-redshift spectrum cannot supply six useful lines under these assumptions. Full-six-line equal-error forecasts therefore remain separate prospective scenarios, not descriptions of this pilot object.

Source code SHA-256: `1e4d20cb3c05dd81e27944eae2664eb43b8a68552f814f071db796dde308faef`.
Source summary SHA-256: `be31314d03f29473e6591ff53a867c7a8163544d3937370e6508081213a42a22`.

The machine-readable checks are in `results/pixel_pilot_validation.json`. The finite-difference and nuisance-projection arithmetic was executed in an independent review script; no source program or observed data were modified.
