# ESPRESSO conventional-model test: analysis scope

Recorded during implementation, before inspecting the fitted differential shifts (2026-09-20). This is an exploratory, publicly model-informed reanalysis, not a blind preregistration.

## Question and hypotheses

In the published right-hand velocity region of the z=1.150793 HE0515−4414 absorber, do six Fe II transitions require different velocity offsets after allowing conventional absorption physics and shared gas structure?

H0: terrestrial isotope wavelengths and abundances with fixed alpha, 45 shared Fe II components, each having column density, velocity and Doppler width; natural damping; Gaussian instrument LSF; pixel integration; per-line continuum normalization and slope; published masks and fixed zero-level corrections.

H1: the same nuisance model plus five independent relative transition shifts (Fe II2374 anchored at zero). It is a generic diagnostic alternative, not a prediction of the Riemann clock or inverse-log-squared law. Common shifts are absorbed in the gas velocities.

## Selection and inputs

Use all six Fe II transitions in the author's right-region fit:2260,2344,2374,2382,2586,2600. Use original f13 core wavelength limits, FITS status==1 and positive errors. No new residual-based clipping. Original author masks and component choices were data informed; explicitly retain that limitation. The original study's free-alpha solution only initializes gas parameters; alpha is fixed at zero in the null forward model. FeII1608 is outside this ESPRESSO spectrum. This is one absorber, not a cosmic-time sample.

Use the pinned ESPRESSO repository commit09ee2cb32fb57ace9cda2943e997930bc7813366, source manifest hashes, published isotope-weighted oscillator strengths without double reweighting, and Gaussian FWHM2.10/2.05km/s in blue/red. Keep the full original flux/error/status arrays. Spectrum is a historical coadd, not a newly obtained2026observation.

## Comparison and controls

Jointly re-fit nuisance parameters under both hypotheses. Check optimizer convergence, bounds and degeneracies. Main weights are the supplied normalized statistical-error array. Repeat with expected-fluctuation weights; repeat under an assumed, explicitly unmeasured AR(1) correlation and an LSF-width perturbation. Check removal of known-blended transitions and numerical subpixel resolution. Compare chi-square improvement against five added parameters; nominal likelihood metrics are conditional on the noise/model and not a universal beyond-conventional significance.

Independent derivative/atomic/mask audit and a controlled injection test must distinguish an insensitive pipeline from a small observed improvement. Never equate absence of detection with proof of immutable constants. No inverse-log-squared fit unless a reproducible residual mode is identified and multiple ages plus physical response make it identifiable.
