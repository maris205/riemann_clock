# Independent UVES cross-instrument implementation and scope audit

Date: 2026-09-21. Reviewer: a separate agent from the adapter/fit author. This review
does not assert a cosmological signal or a complete reproduction of VPFIT.

## Input fidelity

A separate read-only input audit compared the archived original FITS with the
processed NPZ. Wavelength, normalized flux, statistical error, expected-fluctuation
error, continuum and status agree exactly. The validity mask is exactly
`status == 1`, finite flux/error and positive error (279,787 valid pixels in the
full spectrum). The source wavelength coordinate is vacuum heliocentric Angstrom.
The median log-grid spacing is 1.299997181214 km/s; its tiny numerical variation
is floating-point precision.

The fitted 2374/2382/2600 windows contain 156/178/177 pixels, respectively, totaling
511. All are valid. Their wavelength, flux, error and masks agree exactly with
the source samples. The optional transferred ESPRESSO 2586 blend mask changes
none of these three clean-line windows. This confirms mask handling, not the
absence of every unidentified blend.

The archived `CombinedSpecInfo` gives nominal resolving power 53695.928571 in
these bands, correctly rounded to 53696 by the adapter. It also gives arc-line
resolving power 70723.320938. Neither is a measured quasar-illumination line-spread
function (LSF). The nominal fixed Gaussian and fitted per-transition Gaussian
widths are therefore alternative conventional assumptions, not independent LSF
measurements.

The header records 51 exposures and 181652.0812 seconds. Contributing exposure
table dates run from 1999-12-14 to 2009-02-11. Arm/extraction table records must
not be counted as unique exposures.

## Wavelength-correction provenance finding

It would be incorrect to call this coadd simply "uncorrected". The archived UPL
contains 138 arm records, 95 with nonzero `*_VSHT` entries. For example:

```
010_VSHT = 0.278073 0.111905 5800.000000
038_VSHT = 0.439564 0.202713 5800.000000
084_VSHT = 0.139810 0.221430 5800.000000
127_VSHT = 0.307377 0.208591 5800.000000
```

The official [UVES_popler header](https://raw.githubusercontent.com/MTMurphy77/UVES_popler/master/UVES_popler.h)
defines shift in km/s, slope in km/s per 1000 Angstrom, and reference wavelength
in Angstrom. The [wavelength code](https://raw.githubusercontent.com/MTMurphy77/UVES_popler/master/UVES_wpol.c)
applies these fields to the wavelengths. These source links identify the code
semantics; this audit did not rerun every exposure reduction.

Appropriate manuscript wording is: "The public SQUAD coadd records per-exposure
wavelength shift and slope corrections. Exact equivalence to the Kotus et al.
precision-analysis product, and the remaining calibration uncertainty, have not
been established."

## Forward model and numerical checks

The independent executable is `code/cross_instrument_validate.py`. Its current
machine-readable check inventory, source hashes, audited result hashes and exact
numerical discrepancies are in `results/cross_instrument/validation.json`.
Checkpoints from active optimization are not treated as final result products.

The final audit passes **310/310 checks across 17 saved fit products**, including
the three bounded continuation endpoints and both opacity controls. Independent
CGS profiles evaluated at every refined fitted parameter set differ from saved
model flux by at most `4.00122e-10`. Adapter and validator source hashes match
their audited versions. These are implementation checks, not 310 independent
pieces of astrophysical evidence.

The validator verifies:

- Exact null/alternative nuisance parity: 45 shared gas components (135 gas
  parameters), six continuum parameters, three zero corrections and three LSF
  widths; the alternative adds only 2382 and 2600 shifts relative to 2374.
- Exact nested residual equality when those two shifts are zero.
- Pixel-centered quadrature, direct source-array equality and masks.
- Independent frequency-domain Voigt optical depths using CGS constants and
  `scipy.special.voigt_profile`, followed by an independently constructed Gaussian
  kernel, FFT convolution and pixel integration. This differs from the adapter's
  complex-error-function implementation and `gaussian_filter1d` convolution.
- Selected gas, continuum, zero, relative-shift and all LSF-width Jacobian columns
  against centered finite differences, plus mixed directions and AR(1) whitening.
- Recomputed saved chi-square, residuals, Jacobians and model profiles.
- Increased pixel quadrature density at saved fitted parameters.

The optimized LSF derivative capture passes the independent derivative tests.
Its callback changes only the importing process and restores the original
convolution function through `finally`. The appended LSF parameters do not shift
the base model's gas, continuum, zero or relative-shift offsets.

One material numerical issue was found and reported during development: the
first flexible-null result at seven subpixels per UVES pixel changed by
`+1.0508548726` in chi-square when evaluated at 21 subpixels, with squared
noise-weighted model discrepancy `0.0587360858` and maximum pixel discrepancy
`0.0935717207` sigma. That is small in flux but material against a relative-shift
improvement of only a few chi-square units. The adapter author increased the
final-fit sampling to at least 21. Seven-subpixel products are retained as
exploratory records; the validator characterizes their error without labeling
them converged. Refined products receive a further sampling-convergence check.

## Inference limits and required interpretation

The nominal-width fits initially had chi-square near 1400 for 511 pixels and
did not satisfy optimizer convergence. Their approximately positive 80 m/s
relative shifts are therefore poor-model diagnostics, not replicated physical
measurements. Allowing conventional LSF widths reduces chi-square by hundreds
and can change the sign of the fitted relative shifts. This is evidence that
the cross-instrument conclusion depends on resolution modeling.

The final numerical values and optimizer statuses must be taken from the latest
refined products, not those exploratory figures. Saving a finite fit result is
not proof of convergence. Rank deficiency, active bounds and competing minima
must be reported explicitly when present. The diagonal-noise chi-square survival
fields written by the shared base fitter are nominal diagnostics only; they are
not calibrated significance, model-selection probabilities or discovery claims.

The refined statistical-error alternative terminates by `ftol` after 554 function
evaluations, with chi-square 406.8472924 and shifts of -139.2218 and -131.6604 m/s.
Its numerical Jacobian has full rank 149 at the declared SVD threshold. However,
12 parameters are active at bounds: five gas velocities, six gas widths, and
the 2374 LSF width at the lower scale bound 0.65 (equivalent fitted resolving
power about 82609). A full-rank Jacobian therefore does not remove the boundary
and model-dependence limitations. The initial cross-started null has chi-square
418.2560957 after exhausting 400 evaluations; a final bounded continuation
reaches 417.8407498 after another 500 evaluations and remains evaluation-limited.
The original null was 427.1187871 after 700 evaluations. This sequence demonstrates
optimization dependence. The difference from a fully minimized null is not
established.

The selected final endpoint pairs are:

| Error/opacity configuration | Null chi-square | Alternative chi-square | Null termination | Alternative termination |
| --- | ---: | ---: | --- | --- |
| Statistical errors, fixed opacity | 417.8407498 | 406.8472924 | Evaluation limit | `ftol` |
| Expected-fluctuation errors, fixed opacity | 398.1177619 | 388.1353623 | Evaluation limit | `ftol` |
| Expected-fluctuation errors, variable relative opacity | 396.7242284 | 386.5698829 | Evaluation limit | Evaluation limit |

The expected-fluctuation alternative gives shifts -140.2311/-132.2706 m/s;
the provisional opacity alternative gives -139.7512/-130.9429 m/s. All six
selected endpoints have full numerical rank at the declared SVD threshold,
but active parameter bounds and capped null optimization remain. Both fixed-
opacity alternatives and the provisional opacity alternative place the anchor
LSF on its lower width bound. These results support conditional directional
compatibility under the examined models; they do not yield a calibrated
likelihood-ratio significance or an unconstrained physical displacement.

At the refined statistical-error alternative parameters, increasing pixel
sampling from 21 to 49 changes chi-square by +0.00450893, with maximum pixel
change 0.00230932 sigma and squared noise-weighted model change 0.000148245.
This resolves the material seven-subpixel integration issue for this saved
solution. It does not establish optimizer or astrophysical completeness.
The three final continuation products also pass 21-to-49 refinement: their
maximum pixel differences are below 0.00230 sigma and absolute chi-square
changes below 0.00483. Thus the remaining termination/boundary concerns are not
explained by the previous coarse numerical integration.

The source `Notes_FITS_Files.txt` recommends expected-fluctuation errors for model
chi-square fitting. The median expected/statistical ratios in 2374/2382/2600 are
1.08232/1.00710/1.00338. The completed control using the recommended array matters because
it changes the anchor-line weight disproportionately; it still does not provide
the full covariance of the rebinned/coadded data.

Both hypotheses use the same physical model family and nuisance dimensions
except the two shifts, so the comparison is internally nested. It remains
conditional on a 45-component architecture imported from ESPRESSO, unresolved
gas structure, bounded parameter motion, a symmetric Gaussian LSF per line,
unknown residual calibration uncertainty, and noise covariance. Fitted Gaussian
widths do not independently establish the instrument resolution. A fit
improvement from extra relative shifts does not by itself
exclude these conventional explanations.

Finally, UVES and ESPRESSO view the same absorption system. This can test
instrumental persistence of a line-pattern discrepancy; it supplies neither an
independent absorber epoch nor a fit of a cosmic-time law. It cannot test
`1/ln^2(t)` by itself.
