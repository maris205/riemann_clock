# Existing-spectrum analysis: method and interpretation review

Review date: 2026-09-20. This review concerns the two actual UVES SQUAD DR1 spectra already archived in this project. It is independent of the synthetic absorption-line injection pilot. Proposed registration parameters are descriptive summaries of real flux profiles, not measurements of a cosmic Riemann-zero cutoff.

## A defensible empirical contribution now

The data can support a reproducible report of line coverage, continuum noise, observed equivalent widths over specified intervals, saturation/blend warnings, and the stability of **within-absorber relative-profile registration**. The most useful question is whether a simple shared-profile approximation survives changes of line, velocity interval and error model. Failure is an informative design result: it establishes that a photon-noise calculation is inadequate for this dataset. Neither success nor failure alone identifies a new time-dependent physical parameter.

The parent SQUAD release comprises 467 reduced quasar spectra, with varied observing settings and depth. A two-object selection illustrates actual archive conditions but does not estimate their population distribution or supply a controlled redshift trend. Quasar emission redshift and absorber redshift must remain separate. [Murphy et al. 2019](https://arxiv.org/abs/1810.06136)

## Concrete lines and windows

The existing real-flux figure shows that HE 0515−4414 has relatively weak Fe II 1608 and 2374 absorption and deeper 2586 absorption. The 2382 and 2600 troughs approach zero flux. The latter are useful **saturation controls**, not interchangeable versions of the weak lines. The project metadata give the 2374 and 2586 windows the same nominal resolving power; they are the cleaner starting point for a restricted instrumental comparison. The 1608 window has a different nominal resolution, much larger wavelength separation, and lower signal-to-noise, requiring a separate sensitivity test.

The downloaded ±200 km/s windows cover only part of the HE absorber. The published ESPRESSO analysis modelled an approximately 720 km/s complex with 129 components, including very narrow structure. Crucially, that study reports additional unidentified absorption in **Fe II 2586**, also seen in UVES; Fe II 2344 has an identified foreground Ca II blend. Thus a 2374/2586 comparison should be described as blend-sensitive unless the published exclusion regions are transferred and checked. This is a specific source-based reason to avoid claiming a clean physical centroid measurement. [Murphy et al. 2022, sections 3.2, 3.5 and 6](https://arxiv.org/html/2112.05819v1)

The author's Figure 4 and pinned `fit_r.f13` are archived as `data/raw/murphy2022_published_fig4.svg` and `data/raw/murphy2022_fit_r.f13`, with hashes/URLs in `reports/real_analysis_method_sources.json`. The fit input names the broad 2586 wavelength interval 5561.820–5566.100 Å but does **not** specify the manually excised pixels. No exact 2586 blend mask has been invented from this interval or guessed from a plot. The full-window registration remains an uncorrected diagnostic until source masks are explicitly decoded.

For the high-redshift example Q0347−3819, the observed 2586/2600 windows near 1.04 μm have median continuum/error approximately 1–2 and are unsuitable for this precision exercise. The 2344–2382 windows around 0.94–0.96 μm require telluric and saturation screening. The better 1608 window alone does not form a differential same-ion line test. These conclusions follow from the archived flux/error arrays and coverage report; they are not a statement that all high-redshift spectra are inadequate.

Suggested HE sensitivity intervals, all explicitly relative to the same stated reference redshift, are a broad local interval (e.g. −25 to +110 km/s), a central portion (−15 to +40 km/s), and a red portion (+40 to +95 km/s). These are **exploratory analysis choices after inspecting this target**, not preregistered independent experiments. A shift change across these overlapping intervals measures method dependence, not a cosmic-time dependence.

## Wavelength consistency matters already

Recompute every velocity grid from the source vacuum/heliocentric wavelength array and the same NIST Ritz list. Do not use the earlier approximate coverage labels as precision wavelengths. For example, the coverage values 2382.7642 and 2600.1725 Å differ from archived Ritz values 2382.76386 and 2600.17206 Å by about 43 and 51 m/s respectively. The archived Fe II wavelength uncertainties themselves correspond to roughly 13–15 m/s. Shared laboratory calibration/energy-level errors need covariance treatment and do not average down independently with each absorber. [NIST line and uncertainty definitions](https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html)

A common convention is

\[
v=c\ln\!\left[\frac{\lambda_{\rm obs}}{(1+z_{\rm ref})\lambda_{\rm Ritz}}\right].
\]

This is a convenient logarithmic spectral coordinate; identify it as such. Using the first-order wavelength velocity is also acceptable if applied consistently. Even apparently small registration errors should not be interpreted as calibrated absolute velocities.

## What the template fit does and does not mean

A model such as

\[
F_j(v)=C_j(v)\{F_r(v-\Delta v_j)\}^{a_j}
\]

with a local linear continuum, a free depth exponent and a shift is a useful **descriptive registration model**. It approximates the optical-depth scaling expected for one ion. However, observed flux is instrumentally convolved, and in general

\[
\left[L*e^{-\tau}\right]^a\ne L*e^{-a\tau}.
\]

Unresolved saturation therefore creates line-dependent residuals and apparent shifts even with the same nominal resolution and no change of atomic physics. Equal nominal resolving power does not establish identical line-spread functions. An affine-flux alternative provides another model-dependence check, not an independently valid physical model. A physically interpretable velocity shift ultimately requires a latent optical-depth/Voigt model convolved with each instrument response, with blends and alternative velocity structures handled explicitly.

The apparent-optical-depth diagnostic is related but distinct: for continuum-normalized positive flux, `tau_a = −ln F`. Comparison of `tau_a/(f lambda)` for transitions of the same species can reveal unresolved saturation if vetted oscillator strengths and line identifications are available. AOD is not guaranteed to recover the intrinsic optical depth at finite resolution. [Savage & Sembach 1991](https://ntrs.nasa.gov/citations/19910068997) Low-S/N flux noise is also nonlinearly biased by the logarithm, making the poorest near-infrared windows particularly unsuitable for naive AOD treatment. [Fox, Savage & Wakker 2005](https://arxiv.org/abs/astro-ph/0508045)

Specific statistical requirements for the descriptive fit:

- Propagate reference-template noise as well as target noise, including interpolation weights. Interpolation and the original coaddition induce off-diagonal covariance. If a diagonal approximation is used, name it and do not present its error as a complete uncertainty.
- If the propagated variance depends on the fitted exponent/shift, either include the Gaussian log-variance term when claiming a likelihood, or explicitly call the objective a weighted least-squares diagnostic. Avoid interpreting a parameter-dependent weighted sum as a fully normalized likelihood without explanation.
- Report the objective, number of fitted pixels/parameters, residual structure, fit boundaries and all prespecified window/error-model variants. A curvature error measures local numerical precision under the model; it is not a calibrated coverage statement when the profile model is rejected.
- `expected_fluctuation_normalized` and `error_normalized` are distinct source arrays; compare them as a sensitivity analysis. Neither supplies the omitted full covariance.
- Do not silently reject the pixels that prevent a desired shift. Masking from independent identified blends or published exclusion lists is preferable; exploratory residual-based masks must remain labelled exploratory.
- Registration shifts between overlapping intervals/line pairs are correlated. Do not pool them as independent measurements to manufacture a small uncertainty.

Historical UVES/HIRES wavelength distortions are a further physical nuisance, separate from pixel noise; the original supercalibration study found substantial long-range changes. Same-object cross-instrument data are therefore a particularly useful next control. [Whitmore & Murphy 2015](https://arxiv.org/abs/1409.4467)

## What two sightlines cannot establish

There is no defensible new estimate of a universal cutoff index, a cosmic noise floor, or the superiority of an inverse-log-squared law from this pilot. The two absorbers differ in velocity structure, observed wavelength coverage, signal-to-noise and contamination. A redshift difference in raw residuals is confounded by all of these differences. A common atomic response coefficient for each transition has not yet been independently derived for the Riemann hypothesis under study.

The suitable bridge to the existing main paper is therefore: **public laboratory diagnostics constrain what its source experiments show; public astronomical diagnostics determine how a proposed epoch-dependent physical response could be measured.** They are two empirical feasibility layers. Neither layer currently calibrates the proposed cosmic uncertainty amplitude.

## Future design grounded in these actual limitations

1. **Archive stage now:** expand from these examples to a prospectively selected absorber sample with at least three useful same-ion transitions, at least one weak saturation diagnostic, independently defined blends, and traceable exposure calibration. Freeze the inclusion/exclusion criteria before fitting a cosmic response.
2. **Matched present-instrument sample:** prefer absorber redshifts allowing the same line group to be observed. With the previously verified ESPRESSO 380–788 nm range, the six Fe II 1608–2600 lines require approximately `1.3625 <= z_abs <= 2.0306` geometrically. Real throughput, atmospheric bands and detector coverage narrow this range. Reobserve/cross-check HE 0515−4414 as a method target, while remembering that its 1608 transition lies below the ESPRESSO band.
3. **Higher-redshift extension:** obtain the missing long-wavelength weak lines with adequate resolution and continuum signal-to-noise; a broad nominal wavelength range alone will not cure atmospheric contamination. An independent physical transition-response model and wavelength-calibration controls are prerequisites for an age regression.
4. **Conditional ANDES expansion:** its 2026 instrument-team paper describes a 0.4–1.8 μm baseline at resolving power up to 100,000 and a first-light target of 2035. That supports a conditional ten-year extension, not an assured science-survey date or guaranteed velocity precision. [N’Diaye et al. 2026](https://arxiv.org/abs/2607.23555)

For any eventual constant-versus-time-dependent comparison, nuisance constraints and the atomic response must be fixed independently. It is possible for a dataset to constrain a nonzero differential response while still being unable to distinguish inverse-log, inverse-log-squared and nearby slow laws. That distinction should be an explicit separate forecast and held-out test.

## Bibliographic handoff

Verified entries for SQUAD, NIST ASD, the ESPRESSO benchmark, AOD, atomic transition data, historical wavelength distortions and the ANDES 2026 outlook are supplied in `riemann_clock/paper/highz_references.bib`. They are not automatically added to or cited by the main manuscript. NIST version 5.12 is the actual archived snapshot; its publication year is 2024, not the 2026 retrieval date.
