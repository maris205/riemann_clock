# Independent stationarity and local-uncertainty audit

This saved-array audit covers **108 fit records and 1247 checks**, with **0 failures**. It includes recursive `failure_pair` diagnostic records without promoting them to precision points. No nonlinear fit or physical-adequacy claim is added. Per-file hashes and every check are in `results/archive_complete/validation/stationarity_independent.json`.

## First-order definitions

For cost `0.5 * r.T @ r`, let `g = J.T @ r` and use distance to the upper bound for a negative gradient, otherwise distance to the lower bound. `independent_coleman_li` is `max(abs(g) * distance)`. Its unit-box normalization cancels exactly and agrees with the installed SciPy Coleman–Li implementation. The name is appropriate. The separate normalized gradient mapping is coordinate dependent; the column-normalized gradient ignores active constraints and is not a KKT statistic.

Successful recorded stops satisfy `ftol`; they do not establish `gtol=1e-6` or global optimality.

| Canonical selected target | H0 Coleman–Li | H1 Coleman–Li | Delta chi-square |
|---|---:|---:|---:|
| J053007-250329 | 0.1997818623 | 0.6403308834 | 4.17603737 |
| J225719-100104 | 0.2851433755 | 0.2035211875 | 13.53129832 |

## Independent nuisance projection

Independent `gelsd` least squares and pivoted QR both reproduce the implementation. All other gas, continuum, zero-level, line-width and displacement columns are projected out; this is not an uncertainty with those parameters held fixed.

| Target / fit | Canonical | Cutoffs | Rank | Information [(km/s)−2] | Local sigma [m/s] |
|---|---|---|---:|---:|---:|
| J053007-250329 / n16_shift_fromnull_os21 | True | 1e−8, 1e−10, 1e−12 | 66/66 | 59.086060714 | 130.0940643 |
| J225719-100104 / n30_shift_fromnull_os21 | False | 1e−8, 1e−10, 1e−12 | 103/103 | 17.615271799 | 238.26229774 |
| J225719-100104 / n30_shift_recovered_os21 | True | 1e−8, 1e−10, 1e−12 | 103/103 | 19.095537687 | 228.8411145 |

The local nuisance-profiled uncertainty is conditional on the selected architecture and expected-fluctuation diagonal errors. It ignores active-bound truncation, calibration, blends, pixel covariance and the exploratory choice of velocity windows and components. Rank stability does not remove these limitations.

## Sparse profiles and retained history

| Target / point | Canonical | Fixed shift [m/s] | Delta from free chi-square | Fitted / saved parameters |
|---|---|---:|---:|---:|
| J053007-250329 / n16_profile2382_0_os21 | True | 8.59271613 | 3.9065082014 | 66/67 |
| J053007-250329 / n16_profile2382_1_os21 | True | 138.6867804 | 0.97621177358 | 66/67 |
| J053007-250329 / n16_profile2382_2_os21 | True | 268.7808447 | -8.0324130067e-06 | 66/67 |
| J053007-250329 / n16_profile2382_3_os21 | True | 398.874909 | 0.97351285465 | 66/67 |
| J053007-250329 / n16_profile2382_4_os21 | True | 528.9689733 | 3.8896207648 | 66/67 |
| J225719-100104 / n30_profile2382_0_os21 | False | -541.1272318 | 4.1956794649 | 103/104 |
| J225719-100104 / n30_profile2382_1_os21 | False | -302.8649341 | 1.0980127801 | 103/104 |
| J225719-100104 / n30_profile2382_2_os21 | False | -64.60263635 | -2.4514965844e-06 | 103/104 |
| J225719-100104 / n30_profile2382_3_os21 | False | 173.6596614 | 1.0683539535 | 103/104 |
| J225719-100104 / n30_profile2382_4_os21 | False | 411.9219591 | 4.2064416964 | 103/104 |
| J225719-100104 / n30_profile2382_recovered_0_os21 | True | -473.4199343 | 4.3283479639 | 103/104 |
| J225719-100104 / n30_profile2382_recovered_1_os21 | True | -244.5788198 | 1.1018612412 | 103/104 |
| J225719-100104 / n30_profile2382_recovered_2_os21 | True | -15.73770531 | -0.0001399541315 | 103/104 |
| J225719-100104 / n30_profile2382_recovered_3_os21 | True | 213.1034092 | 0.72131966889 | 103/104 |
| J225719-100104 / n30_profile2382_recovered_4_os21 | True | 441.9445237 | 3.7950529727 | 103/104 |

Profile masks exclude only the fixed 2382 coordinate. Pixel arrays, remaining bounds, objectives, degrees of freedom and AICc were independently checked. Any negative difference remains an optimization warning. Five points do not establish global profile minima or interval coverage.

`initial_parameters` was recorded before replacement of the fixed 2382 coordinate. Off-center profile records therefore retain the incoming free-fit seed at that coordinate; the actual active-coordinate starts and final fixed values are correct. No historical records were modified.

## Lower-endpoint optimization caveat

J225719-100104/selected_pair_before_lower_endpoint_recovery.json: incomplete H0 `n30_null_extension_s2_os9` has chi-square 422.187928362, below its selected H0 by 41.53029845 and selected H1 by 32.03767886. Canonical selection: False. These raw comparisons use OS9 versus OS21; exact same-quadrature comparison needs a replay. The original apparent preference cannot be interpreted against the best available conventional endpoint. Retained historical results remain auditable after recovery.

The JSON separates canonical and historical pairs and profiles. Failed-adequacy targets have only stationarity/objective audits here; no displacement uncertainty or cosmic-time observation is created for them.


## Final profile stop status and recovery boundary

The recovered J225719−100104 **+2 local-sigma profile point at 441.9445237 m/s is incomplete** (`optimizer_success=false`, evaluation-budget stop). Its recorded delta chi-square **3.795052973** is the objective at the saved endpoint, not a certified profile minimum. The other four recovered profile points stopped successfully by `ftol`. The center point improves the selected free fit by **0.0001399541** and remains an explicit optimization warning. The profile is asymmetric; it must not be turned into a calibrated confidence interval.

Across all 108 records, **102 stopped with status 2 (`ftol`) and 6 with status 0 (budget)**. Passing numerical-replay checks does not change those optimization statuses. The four newly authorized post-failure diagnostic H0/H1 fits are included in the recursive stationarity checks, but no local uncertainty or precision-point promotion has been performed for them.

For the recovered canonical J225 pair, H0 chi-square is **421.248215119**, H1 is **407.716916802**, and their difference is **13.531298316**. This H0 is below the previously problematic capped null endpoint at 422.187928362. The old higher-basin pair remains in the history tables and the prior optimization warning is not hidden. Finite starts and successful `ftol` termination still do not prove a global optimum.
