# Independent expected-fluctuation sensitivity validation

Passed: **True**; 615 checks; 6 fit products replayed.

The error array was reconstructed directly from original SQUAD FITS primary-array row2, independently of the processed NPZ and sensitivity wrapper. Original row1 statistical errors were also checked. Selected native indices, wavelengths, flux, masks, component counts and parameter labels match the earlier row1 analysis. Residuals and analytic Jacobians are independently reproduced by multiplying the original model outputs by raw row1/raw row2 ratios. All source hashes, saved arrays, objectives, summary arithmetic and selected optimizer termination flags were checked.

|Target|Components|Pixels|Row1 Δχ²|Row2 H0 χ²|Row2 H1 χ²|Row2 Δχ²|Extra shifts|
|---|---:|---:|---:|---:|---:|---:|---:|
|J232128-105122|3|130|6.844288|75.158887|68.998842|6.160046|4|
|J004131-493611|6|138|5.661098|110.035743|104.534085|5.501658|2|

This control changes diagonal pixel weights only. It does not model covariance, reselect gas architecture, calibrate wavelength distortions, resolve all blends or establish a global optimum. Recorded successful stopping is a numerical fact, not proof of physical adequacy. Neither relative shifts nor their χ² improvement directly measure alpha or a cosmological time law.
