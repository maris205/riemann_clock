# Independent validation of the conditional high-redshift forecast

**Result: PASS; 44 numerical checks.** This is validation of the planning arithmetic, not a validation of the physical hypothesis.

The validator does not import the forecast program. It computes cosmic age by numerical integration in scale factor; evaluates the logarithmic law with 60-digit direct powers; projects nuisance modes by least squares; collapses the covariance into six independent bin means; and simulates all six line measurements with explicit offsets, calibration slopes and shared bin errors.

At z = 3, the age is 2.142006656 Gyr and h₂ − 1 = 0.0271017620, or 2.710176% of the assumed present-day amplitude. G₂(3) = 1 by definition. This percentage is not a fractional measured change in any physical constant.

| Conditional scenario | Independent σ(B), m/s | Null 95% half-width, m/s | B for 5σ threshold at 90% power, m/s | Fixed-six-bin floor as N→∞, m/s |
|---|---:|---:|---:|---:|
| UVES_archive_planning | 157.078758 | 307.868708 | 986.698316 | 110.018986 |
| ESPRESSO_narrow_redshift_planning | 110.839167 | 217.240776 | 696.241946 | 68.443041 |
| ANDES_conditional_2035plus | 14.955457 | 29.312157 | 93.943474 | 10.752358 |

The known-covariance GLS is exactly reproduced by six equal-variance bin means: Var(B) = [σ_line²/(m ‖K⊥‖²) + σ_floor²] / Σ_bin(G − mean(G))². The fitted intercept is essential; it removes any constant response offset. The stated 90% power amplitude uses the two-sided Gaussian rejection probability, solved independently by root finding.

The response retains 64.198367% of its squared norm after constant plus log-wavelength slope removal, and 59.169011% after adding curvature. A common shift or pure calibration slope has zero information. If all six bin offsets are free, the cosmic-time column lies in their span and B cannot be identified.

The independent line-level simulation uses 10,000 trials per scenario at zero signal and at the 90%-power amplitude. Large random common shifts and wavelength slopes cancel as expected. All reported intervals and dispersions are statistically compatible with their Gaussian design targets. The complete Monte Carlo results and 99.9% binomial intervals are in `results/forecast_validation.json`.

## Scientific qualifications that must accompany use

- The response K is invented for an experimental-design calculation. It is not an atomic calculation, an alpha sensitivity, or an observed signal.
- The forecast transfers a logarithmic time template to a deterministic line-centroid effect. The paper's extra error standard deviation is a different observable, and no equivalence has been derived.
- The six redshift-bin floors have known Gaussian variances and zero-mean independent draws. Arbitrary free bin offsets or an arbitrary trend-aligned bias are unidentifiable from the proposed signal.
- Equal independent centroid errors, six usable lines per object, equal bin populations, fixed K, and a correctly specified nuisance span are optimistic design assumptions requiring measurement in the pilot.
- The polynomial nuisance is in log wavelength. Real intra-order distortions, blends, isotope mixtures, velocity subcomponents, and model-dependent covariance are not removed by this algebra automatically.
- The Planck-time tstar, fixed exponent, and matter+Lambda age coordinate are assumptions. Normalizing at z=3 removes amplitude dependence, not physical uncertainty about the response law.
- The 5000- or 10000-run Monte Carlo can assess ordinary coverage and approximately 90% power. It cannot establish a 5.7e-7 null false-positive rate; that tail is analytic under the ideal Gaussian model.
- Shape separation thresholds use the expected squared residual against one fixed competitor after refitting its amplitude and intercept. They are not model-selection discovery significances.
- Shape-separation amplitudes of 69 to 1632 km/s extrapolate the toy centroid model far beyond a validated response domain. They are evidence of weak exponent discrimination, not proposed physical target amplitudes.
- The future-instrument scenario is conditional on delivered coverage, sensitivity, suitable objects, and systematic control; its name or numbers do not guarantee a commissioning date or a 2036 result.

## Reproduction

```bash
python code/forecast_sensitivity.py
python code/validate_forecast.py
```

Input summary SHA-256: `544993eb59ebe9bb69216796c4d7e8cf51e35ec65fa2e98e8d3564e00e92f368`.
Forecast code SHA-256: `29b94d6e915f4be33da57c6fa0405e666c85a937f490d54b30f7cbf82c9f622d`.

No observation was fitted in this calculation. A pass supports the internal numerical calculation conditional on its design assumptions only.
