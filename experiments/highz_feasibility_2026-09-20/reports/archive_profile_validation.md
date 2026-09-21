# Independent archive profile implementation validation

Passed: **True**. 3362 checks; 67 saved fits replayed.

Maximum finite-difference Jacobian relative L2 error: 1.12117e-07.

Validation reconstructs native wavelengths, exact indices, errors and validity from the original FITS tar members. It checks isotope oscillator strengths without multiplying abundances twice, all parameter derivatives, an independently evaluated normalized Voigt profile with explicit Gaussian convolution, cache parameter copies, identical H0/H1 pixels, saved objective values, parameter feasibility, and AICc arithmetic.

Quadrature comparisons evaluate the same saved parameters; they do not re-optimize at higher sampling. The complete metrics are in `results/archive_expansion/profile_pilot/independent_validation.json`.

| Fit | χ² at 11 | 33−11 | 55−33 | Largest profile difference / pixel error, 11 vs55 |
|---|---:|---:|---:|---:|
| J004131-493611_n1_null_s2 | 6534.558274 | 0.055559 | 0.004453 | 0.00510268 |
| J004131-493611_n1_shift_s2_fromnull | 4666.240335 | 0.027689 | 0.002226 | 0.00562457 |
| J004131-493611_n2_null_s1 | 502.002410 | 0.000312 | 0.000052 | 0.00784077 |
| J004131-493611_n2_shift_s1_fromnull | 486.586716 | 0.002818 | 0.000253 | 0.00782149 |
| J004131-493611_n3_null_s1 | 341.805873 | 0.000342 | 0.000056 | 0.00853089 |
| J004131-493611_n3_shift_s1_fromnull | 337.064959 | -0.004216 | -0.000308 | 0.00896411 |
| J004131-493611_n4_null_s0 | 197.370940 | 0.000385 | 0.000064 | 0.0103033 |
| J004131-493611_n4_shift_s0_fromnull | 190.688450 | 0.000384 | 0.000064 | 0.0102804 |
| J004131-493611_n6_null_s3_blue_cross | 113.188664 | 0.000386 | 0.000064 | 0.010876 |
| J004131-493611_n6_shift_s3_blue_start | 107.527572 | 0.000387 | 0.000064 | 0.0108632 |
| J232128-105122_n1_null_s1 | 178.760200 | 0.005103 | 0.000429 | 0.00517948 |
| J232128-105122_n1_shift_s1_fromnull | 165.400305 | 0.004946 | 0.000416 | 0.00508391 |
| J232128-105122_n2_null_s0 | 161.755726 | 0.015995 | 0.001300 | 0.00490636 |
| J232128-105122_n2_shift_s0_fromnull | 146.205947 | 0.010740 | 0.000880 | 0.00519265 |
| J232128-105122_n3_null_s1_cross | 80.592853 | 0.000286 | 0.000048 | 0.00589574 |
| J232128-105122_n3_shift_s1_fromnull | 73.748565 | 0.000288 | 0.000048 | 0.00614472 |
| J232128-105122_n4_null_s0 | 76.018730 | 0.000289 | 0.000048 | 0.00589312 |
| J232128-105122_n4_shift_s0_fromnull | 69.412199 | 0.000290 | 0.000048 | 0.00613227 |
| J232128-105122_n6_null_s1 | 66.957703 | 0.002989 | 0.000268 | 0.00583063 |
| J232128-105122_n6_shift_s1_fromnull | 61.457964 | 0.003584 | 0.000312 | 0.00618526 |

## Authoritative selected refined fits

| Target | Components | Δχ² at33 | Δχ² at55, same parameters |
|---|---:|---:|---:|
| J004131-493611 | 6 | 5.661098184 | 5.661098168 |
| J232128-105122 | 3 | 6.844288411 | 6.844288436 |

All selected fits passed their recorded termination status checks. Higher quadrature evaluations are not additional nonlinear fits; source results at33 are the authoritative numerical optimization products.

The pilot remains exploratory: velocity windows and component counts were developed after viewing profiles; diagonal errors are used; several nuisance parameters reach bounds; instrumental and blend uncertainties remain uncalibrated. Successful scipy stopping is not proof of a global optimum. AICc is descriptive under these assumptions, not decisive physical model selection. A parameter called a relative shift is not an alpha measurement or a cosmic time trend.
