# Independent audit: UVES optimizer completion and wider LSF bounds

Date: 2026-09-21. This audit owns only the new validator, this report and
`results/uves_completion/validation.json`. No previous cross-instrument source,
report or result is modified. Initial SHA-256 baselines cover 82 prior files;
the validator checks every baseline before assessing new results.

**Final status:** **209/209 implementation, provenance and saved-result checks
pass**, including the four requested model/bound configurations and one
same-model cross-start. All 82 protected prior files remain byte-for-byte
unchanged. **None of the five fit endpoints meets the predeclared strict
stationarity criteria.** This is a verified limitation of the achieved numerical
solutions, not an implementation-audit failure or a reason to label the
optimization complete.

## Mathematical review

For the weighted least-squares objective

\[
f(x)=\tfrac12 r(x)^T r(x),\qquad g=J^T r,
\]

the bound-constrained first-order conditions are zero gradient in an interior
coordinate, nonnegative gradient at a lower bound, and nonpositive gradient
at an upper bound. This uses the Jacobian of the same weighted residuals used
in the objective.

With finite bounds of positive width \(s=u-l\), define

\[
q=(x-l)/s,\quad g_q=s\odot g,\quad
M=q-\operatorname{clip}(q-g_q,0,1).
\]

For a feasible point, this mapping vanishes exactly at the first-order KKT
conditions. `stationarity()` implements the signs and mapping correctly.
Independent hand-verifiable interior/lower/upper fixtures pass.

The mapping must not be used alone at finite tolerance. For example,
\(q=10^{-12}\), \(g_q=10^{12}\) has mapping \(10^{-12}\), even though a
large gradient points toward its nearby bound. The accompanying Coleman--Li
residual for the bounded trust-region reflective (TRF) method is

\[
C=\max_i\{(x_i-l_i)\max(g_i,0),
                 (u_i-x_i)\max(-g_i,0)\}.
\]

It equals the native final TRF optimality measure in the installed SciPy
implementation. In the example, it is one, so the required joint criteria
correctly reject false stationarity. The validator independently recomputes
this quantity and verifies the saved native metric. If a saved case uses
`dogbox`, its different native free-gradient norm is checked separately.

The fitting script declares stationarity only when native optimality is at
most \(10^{-6}\) **and** the unit-box mapping is at most \(10^{-5}\).
`ftol`/`xtol` success and these stationarity criteria remain separate fields.
Tolerance-projected gradients are supplementary: their active-bound tolerance
can suppress a substantial gradient slightly inside a bound. The validator
also reports normalized primal infeasibility and the parameter with the largest
mapping component.

Passing first-order criteria would establish approximate stationarity for the
specified objective and bounds. It would not establish a global optimum,
exclude saddle points, validate the gas model or LSF, calibrate significance,
or demonstrate a physical change.

## Bounds, nesting and provenance

`CompletionModel.initial()` changes exactly three appended Gaussian-LSF bounds
from scales `[0.65, 1.25]` to `[0.5, 1.4]` in the wider control. Every gas,
continuum, zero correction and relative-shift bound is preserved. The physical
model returns exactly the same residual at a shared parameter vector in both
bound configurations.

For either LSF range, the alternative adds only the two shifts for Fe II 2382
and 2600 relative to 2374. The null and alternative have identical nuisance
parameters and identical residuals when these two shifts are zero. Expected-
fluctuation errors are used throughout this new task. This isolates the wider
LSF range as a specified conventional sensitivity control; it does not make
the allowed range an independently measured instrument response.

The loader's metadata side effect is redirected to the new completion directory
in the fitting script and to a temporary directory in the validator. Old files
remain protected by the pre-recorded hashes. Input snapshots capture the input
record, actual starting vector, bounds, input hash and timestamp. They are
captured provenance, not objects made immutable by filesystem permissions;
the author uses unique case names. The validator checks their saved hashes and
the original input hashes.

## Final numerical checks

`code/uves_completion_validate.py` verifies saved bounds, data/parameter counts,
chi-square, residuals, Jacobians, stationarity classification and source
snapshots. It independently evaluates each fitted profile using a frequency-
domain CGS Voigt calculation and FFT convolution, and compares pixel sampling
21 with 49 at the saved fitted parameters. It launches no optimizer.

The final independent CGS profile calculation differs from saved fitted flux
by at most `4.00842e-10`. Increasing pixel sampling from 21 to 49 changes any
endpoint chi-square by at most `0.00490784`, with maximum single-pixel change
`0.00229079` sigma and maximum squared noise-weighted profile discrepancy
`0.000142857`. Source hashes, starting-record hashes, saved data arrays and
parameter bounds all pass. These checks establish consistency of the saved
calculation, not complete astrophysical modeling.

The parameter with the largest mapping residual at every actual fitted
endpoint also receives a centered finite-difference derivative check. The
analytic columns agree; for example, at `bounded_h0`, the analytic gradient of
`logb_28` is 0.4109253402 and the finite-difference value is 0.4109253357. The
strict nonstationarity is therefore not explained by an erroneous derivative
at the narrow-component solution.

## Verified endpoints and feasible-domain selection

All five runs exhausted their evaluation budgets: 1800 for each primary case,
400 for the cross-start. The independently recomputed diagnostics are:

| Endpoint | Chi-square | Native TRF optimality | Unit-box mapping | Largest mapping parameter |
| --- | ---: | ---: | ---: | --- |
| `bounded_h0` | 397.89334449 | 2.09725102 | 0.96327330 | `logb_28` |
| `bounded_h1` | 388.10859556 | 0.02116666 | 0.04988976 | `logN_27` |
| `wide_h0` | 397.90546138 | 0.03063036 | 0.07158325 | `logN_39` |
| `wide_h1` | 387.77171016 | 0.01579228 | 0.03363581 | `logN_27` |
| `wide_h0_cross_bounded` | 397.85535661 | 0.03194675 | 0.09699201 | `logN_39` |

The cross-started H0 has LSF scales 0.65779677/0.66506429/0.67361508, all inside
the original `[0.65, 1.25]` range. The validator confirms every original bound
is satisfied and that the original model gives **bitwise-identical residuals**
and chi-square 397.8553566088117 at this vector. It is consequently the best
available feasible H0 point for **both** domains. Its improvement over the
separately run bounded H0 must not be attributed to different LSF physics.
No stationarity claim is transferred between domains.

Using that shared feasible H0 point gives descriptive endpoint differences
9.74676105 in the original LSF range and 10.08364645 in the wider range. Neither
is a calibrated likelihood-ratio statistic. The H1 shifts change from
-139.9106/-132.0368 m/s to -144.5814/-136.2845 m/s when widening the range; the
anchor's fitted scale moves from the old lower bound 0.65 to 0.63795825, which
is interior to the wider range. Thus this particular Gaussian-width boundary
does not absorb the relative pattern at the achieved endpoints.

Independent column-normalized Jacobian condition numbers are approximately
`1.53e5–1.60e5` for H0 and `9.50e5–1.08e6` for H1. These quantify weak local
linear parameter combinations despite full numerical rank. They do not prove
global nonidentifiability, and they do not convert an evaluation-limited point
into a stationary point.

The remaining scientific statement is narrowly conditional: the examined
symmetric Gaussian-width range changes the saved relative shifts by several
m/s, while unknown covariance, asymmetric/non-Gaussian LSF, residual wavelength
calibration, gas-component modeling and incomplete minimization remain.
This audit establishes neither exclusion of conventional explanations nor a
cosmic-time dependence.
