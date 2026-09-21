# Independent numerical review: ESPRESSO conventional null

Date: 2026-09-20. Scope: `code/espresso_conventional_null.py`, the published
right-region Fe II component initialization, the actual final ESPRESSO spectrum,
and the implementation of the six-line joint fit. This is an implementation
audit and conditional model comparison, not an independent confirmation of new
physics or an exact reproduction of VPFIT.

The reproducible audit is `code/validate_espresso_null.py`; it writes the
individual assertions, tolerances, measurements, and source hashes to
`results/espresso_null_numerical_validation.json`. The companion mask/covariance
audit is `reports/espresso_mask_covariance_review.md`.

The forward-model audit passes all 64 assertions, including 13 added checks for
the opacity and zero-level nuisance branches. Across all six lines, the
independent physical/convolution calculation differs by at most 4.09e-10 in
normalized flux. The worst selected or mixed-direction Jacobian check has a
relative L2 discrepancy of 9.10e-8. Increasing the subpixel count from 9 to 27
changes the noise-weighted squared model difference by 1.2684e-4, with a maximum
pixel change of 0.00201 sigma; the corresponding chi-square changes by -0.06828
at the perturbed published parameters. Increasing 18 to 27 gives a squared model
difference of 3.10e-6. Quadrature is therefore much smaller than unit-noise
residuals, but objective differences of order 0.07 should not be overinterpreted.

For the added conventional controls, the independent profiles include nonzero
line-opacity factors and nonzero zero-level offsets. Their maximum discrepancy
is 4.05e-10 in normalized flux. Finite differences validate the new opacity and
zero-level derivatives, and the relocated shift derivative columns, to a worst
relative L2 error of 1.76e-8. These checks establish numerical correctness; they
do not establish the adopted nuisance ranges as measured laboratory or
instrumental uncertainties.

## Physical forward model

The optical depth is built before applying the instrumental response. A single
shared set of 45 Fe II component column densities, velocities, and Doppler
widths predicts all six transitions. The isotope entries carry abundance-weighted
oscillator strengths. The normalized flux is `exp(-sum(tau))`, convolved with a
Gaussian instrumental profile and then averaged over subpixels. This ordering
is essential: applying a power to an already broadened template is generally
not equivalent.

The dimensionful opacity coefficient was independently derived from
`sqrt(pi) e^2/(m_e c)` in Gaussian CGS units, with Angstrom and km/s conversions.
It is 1.4973641483561679e-15, agreeing with the implementation to 1.10e-9
fractionally. The damping parameter is `Gamma lambda/(4 pi b)`, including
the corresponding unit conversion. The velocity coordinate is logarithmic,
while the frequency argument retains the exact exponential wavelength ratio.

An independent calculation constructs the absorption cross section in Hz using
`scipy.special.voigt_profile`, physical CGS constants, a separately constructed
Gaussian kernel, and `scipy.signal.fftconvolve`. It also includes the fixed
published zero-level and continuum adjustments. This checks opacity units,
isotope strength weights, frequency-to-velocity conventions, convolution, pixel
integration, and final continuum/zero-level ordering together.

Analytic derivatives with respect to log column density, velocity, log Doppler
width, continuum intercept/slope, and relative transition shifts are checked
against central finite differences. Mixed random directions test parameter
ordering and simultaneous gas/continuum/shift effects. Checks use nonzero shifts
and nonzero continuum perturbations. The independent fine-grid comparison uses
18 and 27 subpixels per data pixel, compared with the fitting default of 9.

## Identifiability and statistical scope

Fixing the Fe II 2374 shift to zero removes the exact gauge freedom between a
common line shift and the shared gas velocities. A transition's shift derivative
equals the sum of the component velocity derivatives in that transition's pixel
block, as required. Any future leave-line-out run must retain this anchor or
explicitly choose a new one. Excluding the anchor while leaving every retained
line shift free recreates the common-shift degeneracy.
The fitting code now rejects the missing-anchor free-shift configuration.

The local gas Jacobian is poorly conditioned. Inverting `J.T @ J` squares its
condition number. The initial implementation used `pinv(J.T @ J)` and discarded
additional weak directions. During review this was replaced by a direct SVD of
`J` with a consistent relative singular-value cutoff of 1e-10. This avoids the
unnecessary normal-equation precision loss. The recorded covariance remains a
local diagnostic: bound-active components, nonlinear poorly identified gas
directions, and alternative gas models still limit the validity of Gaussian
marginal uncertainties and discovery significance.

The parameter-count degrees of freedom and chi-square survival probability are
nominal. Resampled/coadded pixels can be correlated, gas-component selection
comes from a previous fit to these observations, the fit has active parameter
bounds, and continuum/LSF/blend uncertainties are incompletely known. A small
nominal survival probability would not by itself distinguish new physics from
a misspecified conventional model. An acceptable conventional fit provides a
useful negative result: these data do not require extra transition shifts under
the tested conventional model.

Optimizer termination is not proof of a global minimum. Relevant evidence is
the termination reason, residual objective change, scaled optimality, active
bounds, repeated or perturbed starts, and whether the free-shift branch improves
the objective beyond numerical or nuisance-model variations. The alternative
contains the zero-shift null, so a materially worse optimized alternative is an
optimization problem, not physical evidence.

For the final runs the objective tolerance is `ftol=1e-7` and the parameter-step
tolerance is `xtol=1e-8`. A highly degenerate gas model can continue moving its
component parameters long after its predicted flux and chi-square have
stabilized. The defensible convergence criterion for this task therefore
combines paired starts in the nested hypotheses, objective stability, and the
stability of the identifiable relative-shift directions. An `ftol` exit alone
does not certify uniqueness of the component decomposition. Nor does a large
raw or scaled gradient along an almost unconstrained or bounded gas direction
automatically establish a consequential error in the predicted spectrum.

The primary scientific comparison should report both hypotheses' fit quality
and the improvement in chi-square, while retaining the same gas and continuum
flexibility. Restarting the free-shift model from the fitted null and restarting
the zero-shift model from the fitted alternative are useful consistency checks.
Noise-correlation, LSF, continuum, and line-removal variants probe distinct
assumptions. With active bounds and a data-derived component structure,
nominal Wilks/chi-square tail areas are descriptive only; a credible detection
would require calibration under the relevant nuisance and model uncertainty.

This one-absorber analysis does not test a cosmic-time dependence. It cannot
identify the proposed inverse-log-squared law, even if a differential shift
were present. The link between the hypothesized Riemann-resource scale and
atomic transition responses remains an additional, unestablished assumption.
