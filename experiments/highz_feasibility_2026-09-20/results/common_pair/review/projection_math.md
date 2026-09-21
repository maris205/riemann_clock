# Independent common-pair uncertainty mathematics audit

**PASS: 50/50 checks passed.** This audit re-evaluated the two saved expected-fluctuation-row H1 models without running optimization. It makes no edits to fitting code or existing result products.

The requested target is $D = v_{2382}-v_{2374}=1000\,\mathrm{shift}_{2382}$ in m/s, with the 2374 shift fixed to zero. Every other H1 parameter, including all other line shifts, remains a nuisance parameter.

| Target | D (m/s) | Unrestricted local sigma (m/s) | Full / nuisance rank | Normalized full condition number |
|---|---:|---:|---:|---:|
| J232128-105122 | -339.970517064 | 179.582131905 | 33 / 32 | 48.382700 |
| J004131-493611 | 25.635436028 | 327.019043097 | 32 / 31 | 719.892127 |

For the full-J calculation, each Jacobian column was normalized before SVD. Transforming the resulting covariance diagonal back to the original parameter units and multiplying the shift standard deviation by 1000 gives the table. Independently, the derivative with respect to D was projected orthogonally to the space spanned by every other Jacobian column, and sigma was computed as the inverse norm of that projected derivative. Full covariance and nuisance projection agree within $4\times10^{-13}$ m/s. An economic QR projection and a direct inverse of the normalized Gram matrix also agree within $10^{-8}$ m/s.

The saved residuals and Jacobians were replayed from the original physical model, source and coadd hashes matched, and the target-shift Jacobian was checked by finite differences. No chi-square-per-degree-of-freedom rescaling was applied. The errors are conditional on the saved diagonal expected-fluctuation normalization.

Relative SVD cutoffs from 1e-3 through 1e-14 preserve all full and nuisance directions for both targets and leave sigma unchanged. At the deliberately aggressive cutoff 1e-2, J004131 loses two directions: the full-J pseudoinverse reports 315.711428 m/s while nuisance-only projection reports 316.655972 m/s. These truncate different subspaces. The full discarded subspace has a nonzero target-coordinate component, so the finite truncated pseudovariance must not be presented as an unrestricted parameter uncertainty. The full-rank result is 327.019043 m/s.

**Bound and nonlinear caveats.** J232128 has zero_2260 and zero_2374 at their lower bounds; J004131 has zero_2374 at its upper bound. The requested covariance treats these nuisance directions as free in both directions. It therefore describes the unrestricted local tangent approximation to a constrained optimum, not a boundary-aware profile confidence interval. Active-bound directions define a one-sided feasible cone. Artificially freezing them gives a different conditional quantity and does not resolve this inferential issue. The full-rank SVD itself is numerically stable.

The Gauss-Newton Jacobian covariance also omits the residual-weighted second-derivative term in the exact nonlinear Hessian. Neither algebraic agreement nor optimizer termination guarantees global optimality, a parabolic profile, independent noise, adequate gas structure, or calibrated wavelength errors. Formal D/sigma ratios are about 1.89 and 0.078, respectively; these are not validated discovery significances. No alpha variation or cosmic-time law follows.

Machine-readable values, source hashes, cutoff sweeps and checks: `results/common_pair/review/projection_math.json`.
