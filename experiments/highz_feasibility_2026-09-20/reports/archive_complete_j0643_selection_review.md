# J064326−504112 independent pre-fit selection review

Reviewed 2026-09-21 UTC. This bounded audit uses the archived raw FITS, full native-pixel NPZ, metadata, atomic table, previously frozen archive atmospheric policies, and the two supplied selection images. No gas model, displacement, or fitted residual was evaluated.

**Decision: conditionally suitable for an exploratory shared-gas fit to the selected main complex.** Freeze the four lines and velocity window [−145, +100] km/s, retain all source-valid pixels, and state the limitations below. This is not a certification of a blend-free complete absorber or precision wavelength calibration.

All four transitions contain 109 source-valid native pixels with positive finite statistical and expected-fluctuation uncertainties: 436 / 436 total. All statuses equal 1. Direct extraction from the raw FITS reproduces flux, both uncertainty arrays, continuum, status, and the full logarithmic wavelength grid; the raw FITS SHA-256 matches archived metadata. The earlier line-window NPZ files select the same indices. Native spacing is 2.24999156 km/s.

| Fe II | Observed retained wavelength (Å) | Median CNR, expected error | Min flux | Pixels F < 0.1 | First retained flux ± expected error |
|---|---|---:|---:|---:|---|
| 1608 | 5882.480–5887.250 | 56.2 | -0.006 | 17 | 0.9889 ± 0.0179 |
| 1611 | 5892.555–5897.333 | 50.1 | 0.810 | 0 | 0.9895 ± 0.0193 |
| 2374 | 8683.994–8691.035 | 23.9 | -0.032 | 15 | 1.0211 ± 0.0527 |
| 2382 | 8714.353–8721.419 | 37.0 | -0.034 | 39 | 0.9713 ± 0.0379 |

## Wavelength placement and source quality

The quasar emission redshift in the FITS header is 3.09. Its Lyα wavelength plus the archive 3000 km/s buffer is 5022.096 Å. All four windows lie redward of this and outside the archived conservative telluric and broad-water bands, including ±30 km/s heliocentric padding. This rules out those specific broad-band flags, not weak atmospheric or unrelated absorption.

Fe II 1611 retains a strong-sky proxy warning: its 5892.555–5897.333 Å actual window approaches/overlaps the conservative 5891 ± 1 and 5897 ± 1 Å sky ranges when padded. Fe II 1608 no longer overlaps those ranges after narrowing to the present window. The original ±250 km/s readiness windows flagged both lines. There is no exposure-specific sky diagnosis in this audit. The 1611 warning must be retained in subsequent reporting, and a with/without-1611 sensitivity fit is recommended.

## Boundary and absorption assessment

The actual first retained pixels of all four lines are near unity; the −145 km/s boundary is a defensible local separation point. It is not justified to say that this window covers an isolated complete absorption system. A 2382 feature around −160 km/s lies just outside the blue boundary, and only a few native pixels separate it from the main onset. The 1611 spectrum has deeper features around −170 and +135 km/s outside the chosen interval. Their identities have not been established: neither call them proven blends nor assume they share the selected gas model. The 1608 surroundings are depressed by roughly 5% in places, which also leaves continuum and weak-absorption ambiguity.

The mean flux over [−145, −137] km/s is 0.9508 ± 0.0089 for 1608 and 0.8046 ± 0.0177 for 2382, using diagonal expected errors. These low slice means largely include the onset of absorption already inside the window; they do not by themselves demonstrate that the first boundary pixel cuts strong absorption. The first-pixel values above and the supplied images distinguish those cases. The red boundary near +100 km/s returns close to the continuum, but nearby weak structure remains. A boundary sensitivity check is appropriate if a reported relative shift depends on the edge.

## Value of Fe II 1611

The weak 1611 line is suitable as a saturation constraint: its minimum retained flux is 0.810, with no pixel below 0.1, while the other three transitions have deeply saturated regions. In the common central [−45, +40] km/s interval, its velocity equivalent width is 7.015 km/s with a diagonal expected-error signal-to-noise ratio of 24.66. This is a flux-detection diagnostic, not a physical-variation significance. Its broad weak absorption overlaps the central complex, making its inclusion scientifically useful despite the sky-feature caveat. The full fixed interval gives 10.094 km/s and a diagonal signal-to-noise ratio of 20.92.

Use the isotope-weighted laboratory references and summed isotope oscillator strengths already supplied by the atomic table; its individual isotope strengths are abundance weighted. The weak line must constrain the same nonlinear opacity model, rather than supply a separate arbitrary absorption profile.

## Limits and reproducibility

The catalogue absorption redshift 2.659 is rounded and defines a plotting coordinate, not a precision systemic redshift. The selection was made after inspecting profiles and is exploratory. Pixel covariance, per-exposure calibration, detailed sky/telluric response, line-spread function uncertainty, isotope-mixture uncertainty, and identification of nearby features are not resolved by this audit. No alpha or cosmic-time relation has been fitted.

Machine-readable references, wavelength positions, edge statistics, raw-source checks and SHA-256 hashes are in `results/archive_complete/J064326-504112/selection_independent_review.json`.
