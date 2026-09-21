# Independent review of observed-profile metrics

Status: PASS for the descriptive calculation and stated limited use. This is not an astrophysical transition-shift or cosmic-time inference.

This independent review read both full-spectrum NPZ files and the archived NIST selected-transitions CSV directly. It did not import or execute functions from `measure_observed_profiles.py`. The numerical comparison used a separately implemented scalar optical-velocity integration and checked all 12 result rows.

## Definitions and derivatives

For each pixel, construct edges from neighboring full-source velocity centres, clip the interval to the common target window, and denote retained width by d_i and retained midpoint by m_i. Using normalized flux F_i, define A = sum[d_i(1−F_i)], B = sum[d_i m_i(1−F_i)], rest EW = λ_rest A/c, and profile centroid C = B/A. Partial edge pixels require the midpoint of the retained interval; no positive clipping of the deficit is permissible. These formulas integrate piecewise-constant pixels exactly in the explicitly defined optical-velocity coordinate.

The gradients with respect to measured normalized flux are ∂EW/∂F_i = −λ_rest d_i/c and ∂C/∂F_i = d_i(C−m_i)/A. Conditional independent-pixel variances are the sums of squared gradient-times-error products. A covariance model must instead use gᵀCov g. Full covariance is not delivered by the archived combined-spectrum product.

A continuum multiplier 1+s means F_i→F_i/(1+s), not addition of an arbitrary absorption offset. With full coverage the first derivative of centroid with respect to s is Δv(v_window_mid−C)/A. Incomplete windows require the actual retained-pixel sums. The script correctly applies the finite perturbations and keeps these sensitivities separate from statistical confidence intervals.

AR(1) stress tests use original source-pixel indices, retaining any gaps, with rho=0.3 and 0.6. These are assumed covariances, not covariance estimates. The production script also reports the archive's expected-fluctuation array as a separate diagonal calculation; this array is not a replacement for a full covariance matrix.

## Independent numerical comparison

Checked all deterministic output fields relating to EW, centroid, diagonal errors, expected-fluctuation errors, continuum perturbations, four range perturbations, AR(1) sensitivities, coverage, pixel count, and flux diagnostics. All agree within absolute 1e-10 in their stored units and relative 1e-11. All 12 windows have full archive-valid coverage to floating-point precision; archive validity is not a telluric or blend clearance.

| Target | Transition | EW mÅ | Conditional EW sigma mÅ | Independent centroid km/s | Conditional centroid sigma km/s | Use |
|---|---|---:|---:|---:|---:|---|
| J051707-441055 | Fe II 1608 | 179.198016 | 0.499029 | 33.883125 | 0.117074 | descriptor |
| J051707-441055 | Fe II 2344 | 506.050987 | 0.345683 | 35.305623 | 0.030842 | descriptor |
| J051707-441055 | Fe II 2374 | 220.838471 | 0.435609 | 33.871420 | 0.082706 | descriptor |
| J051707-441055 | Fe II 2382 | 757.107239 | 0.263795 | 35.606734 | 0.017769 | descriptor |
| J051707-441055 | Fe II 2586 | 442.981497 | 0.377560 | 34.623453 | 0.037435 | descriptor |
| J051707-441055 | Fe II 2600 | 777.266474 | 0.308179 | 35.346136 | 0.019735 | descriptor |
| J034943-381030 | Fe II 1608 | 235.568127 | 0.651138 | -20.194088 | 0.108141 | descriptor |
| J034943-381030 | Fe II 2344 | 585.765475 | 1.807344 | -25.733456 | 0.117824 | descriptor |
| J034943-381030 | Fe II 2374 | 353.928007 | 5.001242 | -17.600910 | 0.526213 | descriptor |
| J034943-381030 | Fe II 2382 | 645.035927 | 3.719600 | -16.927191 | 0.220443 | descriptor |
| J034943-381030 | Fe II 2586 | 477.298969 | 76.417793 | -36.279815 | 6.304631 | excluded; independent diagnostic only |
| J034943-381030 | Fe II 2600 | 332.613439 | 47.799347 | 0.526853 | 6.786456 | excluded; independent diagnostic only |

The final two high-redshift centroid numbers above are shown only to document independent arithmetic. Their production CSV/JSON centroid fields are correctly absent/null, and they must remain excluded from scientific centroid tables, figures and inference.

An additional independent Monte Carlo check used seed 884132 and 20,000 independent Gaussian draws for each of the ten non-excluded profiles, about their observed flux values. The empirical-to-linearized centroid-sigma ratios ranged from 0.99031 to 1.00592. Symmetric finite differences of every flux derivative (step 1e-6) agreed with the analytic centroid gradient to less than 9.4e-9 km/s per normalized-flux unit. These checks establish arithmetic and first-order noise propagation in the stated conditional model; they do not validate the unmeasured covariance, calibration or line formation model.

## Scientific interpretation and required cautions

- The integration windows were chosen after inspection and are exploratory, common within each target. EW values represent the selected portion of an absorption complex, not necessarily its total EW.
- The reported velocity zero is set by an approximate reference absorber redshift. Whole-profile centroids are not absolute atomic rest-frequency measurements.
- Saturation, velocity-component weights, blending and instrumental effects can produce different centroids in different transitions. Their individual contributions have not been decomposed. The raw differences do not justify a common-centroid chi-square test of new physics, a varying-alpha fit, or a cosmic-time fit.
- The high-redshift 2344/2374/2382 Å windows remain screening quantities with atmospheric-water contamination risk. No telluric correction has been inferred from their archive-valid masks. The 2586/2600 Å windows remain excluded for centroid work.
- The laboratory rest-wavelength uncertainties, their unknown cross-line correlations, wavelength calibration and astrophysical profile uncertainty are not included in the pixel-only error bars. Laboratory values recorded in the result are therefore contextual information, not an already combined uncertainty budget.
- Continuum perturbations and window choices can move descriptive centroids more than the formal pixel errors. This is a useful existing-data result and a reason to require component modelling and calibration before a physical shift analysis.
- The measured wavelengths, fluxes and uncertainties are processed, continuum-normalized archive products, not newly acquired raw photon measurements.

## Reviewed file identities

- `code/measure_observed_profiles.py`: SHA-256 `3126217b97f9686ac0ff9a1ebb0849ad62d5770858de4370fdd6426d6b826860`
- `results/observed_profile_metrics.csv`: SHA-256 `71e66d6315b8acd5a58fa71f081c2cb4002d7108c34dceadab765a1f4a0c244b`
- `data/atomic/selected_transitions.csv`: SHA-256 `57fd0c08d6e4b59922443c0b14d4014ee439193fc57403ba099c2f7da474e02f`
- `data/processed/J051707-441055_squad_dr1.npz`: SHA-256 `8d70330f95090cea5955c64b662a9955b5de0cab95b031ba79582263a604b194`
- `data/processed/J034943-381030_squad_dr1.npz`: SHA-256 `007f300f0c803c7ff5e5c309bf76338ffe3e15f09b94fe746a41e28973c22c5d`

- `results/observed_profile_metrics.json`: SHA-256 `c1e248bfdd17eb6b507e6dd6d203f1212cd23636dfc24fbb8d3f7bb7439110a4`

Max absolute discrepancies by field (stored units):

- `centroid_continuum_1pct_max_change_km_s`: 0
- `centroid_error_assumed_AR1_rho0p3_km_s`: 6.94e-18
- `centroid_error_assumed_AR1_rho0p6_km_s`: 5.55e-17
- `centroid_error_diagonal_km_s`: 0
- `centroid_error_expected_fluctuation_diagonal_km_s`: 0
- `centroid_km_s`: 0
- `centroid_range_10kms_max_change_km_s`: 0
- `coverage_fraction`: 0
- `flux_min`: 0
- `flux_percentile05`: 0
- `fraction_pixels_flux_below_0p1`: 0
- `fraction_pixels_flux_below_3sigma`: 0
- `median_continuum_to_error`: 0
- `rest_EW_continuum_1pct_max_change_mA`: 0
- `rest_EW_error_assumed_AR1_rho0p3_mA`: 4.44e-16
- `rest_EW_error_assumed_AR1_rho0p6_mA`: 1.11e-16
- `rest_EW_error_diagonal_mA`: 0
- `rest_EW_error_expected_fluctuation_diagonal_mA`: 0
- `rest_EW_mA`: 0
- `rest_EW_range_10kms_max_change_mA`: 0

Final metadata refresh: hashes above cover the final regenerated script and CSV/JSON products. The additional expected-fluctuation validity flag and plot bounds do not alter the reviewed numerical measurements. The report now states that the observed profile differences can be produced by the listed astrophysical/instrumental effects without asserting that their individual causes were identified; this wording concern is resolved.
