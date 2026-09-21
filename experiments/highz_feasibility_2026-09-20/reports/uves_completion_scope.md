# Bounded completion of the UVES conditional diagnostic

This is a follow-up informed by the previous exploratory result, not a blinded
preregistration. Previous cross-instrument files are preserved. The reason for
this follow-up is specific: the reference-line Gaussian width reached its allowed
lower boundary, and the conventional null endpoints reached evaluation limits.

## Fixed question and data

The three existing UVES Fe II windows (2374, 2382, 2600; 511 pixels) and their
archive-recommended expected-fluctuation errors remain the data. The published
45-component architecture and terrestrial isotope atomic table remain fixed as
model inputs, while all gas/continuum/zero/width nuisance values are refitted.
H0 fixes laboratory relative transition positions; H1 adds only two relative
velocities with 2374 as anchor. No alpha or cosmic-time parameter is inferred.

## Four primary cases

1. H0 with the original Gaussian FWHM multiplier interval [0.65, 1.25].
2. H1 with the same interval and the same nuisance freedom.
3. H0 with the wider interval [0.5, 1.4].
4. H1 with that same wider interval and nuisance freedom.

The wider interval is a sensitivity range, not an independently measured LSF
prior. All other bounds and pixel selection are unchanged. Each primary case
starts from a captured complete record of the prior selected endpoint. Unique
case names prevent overwriting these snapshots in this run. Numerical source
hashes are captured at module import, together with exact library versions.

Use 21 subpixels per original pixel, with a final independent denser-grid check
at each selected endpoint. The strict solver tolerances are ftol=1e-10,
xtol=1e-11, gtol=1e-6; the initial bounded budget is 1800 function evaluations per
case. A cross-start from an opposing-hypothesis endpoint is allowed if required
to reveal a lower null basin. Do not interpret budget exhaustion as convergence.

## Termination and interpretation

Record the optimizer message and status separately from first-order diagnostics.
For f = (1/2)||r||^2 and g = J^T r, define q=(p-l)/(u-l) and gq=(u-l)g. The
independent unit-box gradient mapping is ||q-clip(q-gq,0,1)||_infinity. Report raw
feasibility, active bounds, the native TRF Coleman-Li scaled gradient norm, this
mapping, and projected derivatives with correct KKT bound signs. The declared
strict stationarity flag requires both native optimality <=1e-6 and unit-box
mapping <=1e-5. It is a local criterion, not a certificate of a global minimum.

If no stable stationary comparison is obtained within this bounded task, retain
that result explicitly and characterize the basin/LSF dependence; do not turn a
finite endpoint difference into a calibrated significance. Shared laboratory
inputs and architecture remain common between UVES and ESPRESSO despite
independent photons. The same absorber does not create another cosmic epoch.

Asymmetric instrumental shapes would require a separately justified calibration
model. They are not part of these four matched cases, and no arbitrary new
shape grid is planned. The data, covariance, residual wavelength calibration,
and gas-model uncertainties limit any beyond-conventional conclusion.

## Recorded cross-start decision

After `bounded_h0` reached its 1800-evaluation cap at chi2=397.8933444935101,
the concurrently continued wider H0 had a higher objective, although that
bounded-H0 vector is feasible in the wider domain. One additional same-model
start, `wide_h0_cross_bounded`, therefore uses the captured `bounded_h0.json`
endpoint with the wider bounds and at most 400 further evaluations. This
implements the allowed cross-start to distinguish optimization history from
the effect of widening a feasible domain; it introduces no physical nuisance
or alternative model. Its stop and stationarity flags are reported separately.
