# Independent validation of adjacent-order residual diagnostics

PASS: 14533 / 14533 checks passed; 0 failed.

The review reconstructs all same-transition paired order shifts and variances directly from the archived fit covariance. The common Fe II 2374 anchor cancels exactly and contributes no extra variance. Lower and upper order native rows, and the two traces, are disjoint.

All 136 new one-row fits were replayed with independent wavelength quadrature-coordinate assembly, QR continuum elimination and centered derivative covariance. Independent bounded scalar optimizations cover every row for exposures 0, 8 and 16, plus every boundary fit. Source hashes include all frozen inputs and all 17 raw FITS files.

Point metadata, continuous inverse-interpolated reference positions and fractional pixel phases were reconstructed. The phase is not the remainder of an integer pixel index. Fixed-template strength categories, velocity bins and phase-group residual moments were rebuilt from model minus native flux divided by native ERRDATA.

All weighted means, heterogeneity statistics, epoch contrasts, seven predictor regressions, four paired sine/cosine phase regressions, all leave-one-out means/slopes, predictor correlations and trace contrasts were independently reconstructed. No residual chi-square variance rescaling is used.

Maximum relative centered-derivative covariance difference: 1.13163e-06.
Maximum independent scalar-refit shift difference: 0.00102134 m/s.

One exposure-10 Fe II 2382 lower-order trace-0 fit reaches the +1.5 km/s velocity limit and remains in the analysis. Its symmetric local covariance is not a calibrated interval; the associated 2382 trace statistics must be treated descriptively.

The additional 68 common-original-support fits were also replayed with independent native-coordinate assembly and nuisance elimination. Every retained mask index and opposite-order original interval union was reconstructed by a separate connected-component algorithm; all original source hashes remain unchanged. The geometric control retains 62,154 of 62,404 pixels (2382: 30,995 / 31,182; 2600: 31,159 / 31,222). All paired means, covariance sums and descriptive changes reproduce. Final native-bin supports remain slightly unequal near edges/holes, and no independent new-minus-old error is assigned to overlapping estimators.

Scientific limits: the protocol was frozen after the original approximately −60 m/s result and is exploratory. BERV, observing epoch and detector coordinate are strongly confounded. SNR derives partly from the fitted data and is a descriptive proxy. Fixed reference phase does not represent every blended component. All native error and model assumptions remain conditional; matching traces or exposures does not independently confirm physical drift. No cosmic-time or discovery claim follows from this audit.

Detailed checks and exact file hashes: `results/order_response/residual_validation.json`.

Failures: []
