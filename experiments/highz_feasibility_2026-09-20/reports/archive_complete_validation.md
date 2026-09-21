# Independent bounded archive campaign validation

Passed: **True**. 6469 checks; 108 stored fits; targets: J053007-250329, J064326-504112, J225719-100104, J233156-090802.

No final fits yet checked for: none.

The audit reconstructs original native wavelength indices and all likelihood arrays from archived FITS, explicitly rejects repeated global native pixels across likelihood regions, verifies frozen H0/H1 pixels and bounds, finite-differences every parameter family, and computes representative flux profiles independently with a cgs optical-depth coefficient, scipy Voigt function and explicit Gaussian convolution. Saved objectives/Jacobians/arrays, stopping flags, AICc selection and conditional adequacy gates are replayed.

Maximum finite-difference relative L2 derivative error: 2.168030981404843e-08.

|Target|Selected fit|χ² at21|χ²41−21|Maximum flux difference / error|
|---|---|---:|---:|---:|
|J053007-250329|n16_null_refined_os21|482.298466|1.64009e-05|0.00103453|
|J053007-250329|n16_shift_fromnull_os21|478.122428|1.52676e-05|0.0010316|
|J064326-504112|n16_null_refined_os21|532.225274|8.30672e-06|0.00130737|
|J064326-504112|failure_n16_shift_from_primary_os21|527.587723|8.1276e-06|0.00131394|
|J225719-100104|n30_null_recovered_cross_os21|421.248215|1.02007e-05|0.000704168|
|J225719-100104|n30_shift_recovered_os21|407.716917|1.03161e-05|0.000667776|
|J233156-090802|n26_null_refined_os21|1672.841768|1.79965e-06|0.000217958|
|J233156-090802|diagnostic_H1_from_primary_os21|1667.298487|1.42292e-06|0.000218272|

Quadrature checks reuse fixed fitted parameters; they are not new optimizations. Successful ftol/xtol stops do not imply gtol stationarity or a global optimum. The separate stationarity/projection audit records first-order values and local conditional nuisance-projected errors. These are not calibrated discovery significances or alpha/time-law measurements.

Sparse-profile provenance: the original generic driver saved `initial_parameters` before inserting the fixed2382 value. The actual full initial vector is reconstructed by setting that one entry to `fixed_2382_shift_km_s`; active-coordinate starts and fit results are unchanged. This audit preserves source/results and records the caveat.

The actual per-count selector uses the lowest chi-square successful stop when available, otherwise the lowest incomplete endpoint. The latter fallback is recorded explicitly in `validation/per_count_selection_audit.json`; every lower unfinished endpoint is retained. Later bounded recovery fits are audited as additional products, without rewriting original selection history.

|Target|Stage|Count|Selected endpoint|Lowest unfinished endpoint|Selected minus unfinished χ²|
|---|---|---:|---|---|---:|
|J225719-100104|bounded_extension_count_grid|30|n30_null_extension_s1_os9|n30_null_extension_s2_os9|41.530041|
