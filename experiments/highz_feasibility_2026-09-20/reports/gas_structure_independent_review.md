# Independent audit of merged gas-structure sensitivity

Generated UTC: 2026-09-20T17:34:05.471607+00:00

Numerical result: **333/333 PASS**.

The adapter uses fixed groups of the published 45 centroids. The 1.025 km/s adjacency rule produces 40 components. Independent scalar sums reproduce total column, column-weighted velocity and the finite Doppler-core second moment; no finite Lorentzian/Voigt moment is asserted. Bounds remain centered on this published merged initializer.

Tests include exact expansion of every merged cloud into coincident subclouds in the unchanged 45-component engine, equality of profiles and folded analytic Jacobians, the no-merging limit, scalar optical-depth and explicit Gaussian convolution calculations, finite differences of gas/continuum/shift columns, and direct covariance inversion using true pixel separations. An adversarial input reverses fitted centroid order to verify that initialization does not regroup fitted centroids.

4 saved optimizer endpoints replayed. Final comparison present: True. Numerical replay is distinct from optimizer convergence. All retained statuses and stopping messages are listed below.

| Endpoint | chi2 | Success | nfev | Optimality |
|---|---:|---|---:|---:|
| merged40_alternative_cross | 2249.24125116 | True | 145 | 4.01534 |
| merged40_alternative_from_null | 2236.29029467 | True | 27 | 2.19628 |
| merged40_null_cross | 2255.15580438 | True | 39 | 5.17792 |
| merged40_null_published | 2255.20363926 | True | 146 | 6.27398 |

Scientific limits: this is an exploratory sensitivity protocol frozen before its new refits, not a prospectively preregistered discovery test. One predetermined component-merging architecture cannot exhaust gas-structure uncertainty. The changed architecture also changes centroid bounds through its initializer. A formally successful local least-squares termination is not proof of a global minimum. Empirical continuum covariance transfer and wavelength calibration remain conditional. Comparisons with the 45-component model are a combined structure/initialization sensitivity, not a nested significance test or evidence for a cosmic-time law.

Full checks: `results/gas_structure/independent_validation.json`.
