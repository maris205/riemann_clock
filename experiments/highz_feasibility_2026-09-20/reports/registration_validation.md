# Independent validation of observed-profile registration

**PASS: 156 checks of the final 24 fits, plus 15 checks of the data-informed design.** Actual observed flux is used, but this remains a descriptive registration analysis and not a measurement of physical constant variation.

The review independently read the three archived windows and NIST wavelengths, reconstructed velocity coordinates and noisy-template predictions, propagated template errors with a full interpolation matrix, and reproduced all stored objectives and residual statistics. The analysis program was not imported. A third random seed in a broad independent search checked every final stored minimum; no materially lower objective was found. These numerical searches do not certify global optimality. Local curvature was separately computed from analytic-gradient derivatives with a smaller step for all retained formal errors.

| Case | Registration shift, km/s | Reported local σ, km/s | Residual statistic / nominal dof |
|---|---:|---:|---:|
| 2586/whole/fluctuation/power | 0.016660 | 0.012269 | 5.262911 |
| 2586/central/fluctuation/power | 0.022540 | 0.013915 | 7.112262 |
| 2586/red/fluctuation/power | -0.034293 | 0.027389 | 0.880368 |
| 2600/whole/fluctuation/power | -0.180583 | withheld | 34.172412 |
| 2600/central/fluctuation/power | -0.180548 | withheld | 49.109191 |
| 2600/red/fluctuation/power | -0.547266 | withheld | 10.788468 |

Initial local-minimum problems in the 2600 red power fit and four central/red affine fits were corrected by running two broader searches for all 24 cases. The final results passed a third independent-seed search. Interpolation knots had also produced spuriously small Hessian errors; the final source now withholds errors at knots, parameter bounds or optimizer failure. The primary 2586 whole/red conclusions remain unchanged.

## Interpretation and limits

- The sign is correct: m(v)=T(v−δ) places a template feature at v_feature+δ, so δ>0 points towards larger/redder velocity.
- The interpolation variance and transformed-template derivative are correct to first order in template noise under independent input pixels. This is not an exact marginal likelihood for a latent absorption profile.
- The objective includes the parameter-dependent log variance term. Its Hessian is a minus-twice-log-likelihood Hessian, so the local covariance convention is 2 H⁻¹.
- Linear interpolation creates cross-pixel covariance even for independent template inputs. The adopted objective retains only its diagonal. For whole 2586 / fluctuation / power, the omitted adjacent correlations have median about 0.289 and maximum about 0.317 after adding target noise.
- Convolved-flux powers do not commute with line-spread-function convolution. The published 2586 blend, 2600 saturation, continuum choices and fit region can move registration parameters without physical frequency drift.
- The primary whole 2586 residual statistic per nominal degree of freedom is about 5.263, while the selected red subregion is about 0.880. The latter was selected after inspecting the data and is neither independent confirmation nor a cosmological test.
- The interpolation-knot, boundary and optimizer guards correctly withhold unstable curvature errors. These guards do not repair model inadequacy or make retained local formal errors calibrated confidence intervals.
- The residual statistics use nominal degrees of freedom. Their ratios are diagnostics, not calibrated goodness-of-fit p values, because correlations, noisy-template effects and model discrepancies remain.
- Laboratory/calibration covariance, transition-dependent instrumental profiles and the extended absorber structure are omitted. No reported registration parameter is an alpha, cosmic-aging or Riemann-cutoff estimate.

Source code SHA-256: `9ad4120aa7a5055f2fa6c9413708c8c8a34b839363c566e720480dde8f82f515`.
Source results SHA-256: `6cee0542125c5daecd79710d462d4a1fe68c9011fc071277bc54cf44d956d47f`.

Per-case checks, independent searches and interpolation-correlation diagnostics are recorded in `results/registration_validation.json`. No astronomical detection significance is assigned by this review.

## Supplement: data-informed design

The rank/null demonstration and exposure ratios in `design_from_observed_data.py` and `results/data_informed_design.json` passed 15 additional checks. The nuisance matrix has rank 4, certified by a nonzero 4 × 4 minor. The signal obeys the exact relation

`signal = G_low × K_column + (G_high − G_low) × K_1608 × high_z_redshift_column`.

Thus the full 7-row, 5-column design still has rank 4. Even granting six usable low-redshift lines, a single high-redshift line cannot identify B when its absorber redshift and the response-mode intercept are free. This is a structural degeneracy; more exposure on that one transition cannot remove it.

The assumed time-law leverage is ΔG₂ = 0.494575937; arbitrary illustrative amplitudes B = 100 or 30 m/s correspond to mode differences 49.457594 or 14.837278 m/s. These are conditional algebraic examples, not measured or predicted physical effects.

Every exposure factor independently reproduces `max(1, (target CNR / archived CNR)^2)`. For the two excluded high-redshift windows, the numerical factors are:

| Transition | Current archived CNR | Factor to CNR 30 | Factor to CNR 50 |
|---|---:|---:|---:|
| FeII_2586.6493 | 1.311821 | 522.990 | 1452.749 |
| FeII_2600.1725 | 2.063353 | 211.395 | 587.209 |

These factors assume unchanged observing conditions and photon-noise scaling. They do not specify telescope hours, account for overheads, or cure atmospheric absorption or blends. Nominal wavelengths in this selection table serve coverage planning, while the actual-profile analysis separately recalculates precise NIST coordinates.

Design code SHA-256: `4b3cc70d2236e3250c7dd94d197fab267d632c21eadd6f22614a88fdc620b391`.
Design results SHA-256: `81990b976b395e95cd79c4ed32e67911cc14ffb94704c823504ac0c4dfea46fa`.
