# Independent static review of common-pair refinement

**Conclusion:** No blocking source-logic error was found for the currently saved run. The full numerical validator remains responsible for replaying every endpoint, masks, derivatives, optimization status and file integrity. This review did not run optimization or modify fitting code. It examined both companion source files and read the completed refined summary.

## Profile construction and nuisance treatment

`common_pair_fit.fit` fixes only `shift_2382` when a profile coordinate is supplied. The full vector is reconstructed from the fixed coordinate and the other optimized coordinates, and the Jacobian is restricted to exactly the optimized columns. Other line shifts remain nuisance parameters; 2374 is the existing fixed anchor. H0 instead instantiates `ExpectedModel(..., free=False)`, eliminating every relative line shift while retaining the frozen gas architecture, windows, line list, expected-fluctuation errors and all conventional nuisance types. Thus the D=0 profile and full H0 are correctly distinguished.

Each refined grid coordinate receives a warm-start fit and an H0-derived start, retaining the lowest objective among all saved starts at that same coordinate. Selection deliberately includes unfinished endpoints and preserves optimizer status; this avoids hiding a lower unfinished endpoint but requires the downstream status check. The completed summaries report that every selected profile point terminated. A finite collection of starts does not certify the global profile minimum.

The companion selects the lowest recorded endpoint, then explicitly releases the fixed coordinate when its objective is below the best free endpoint by more than 1e-5. The final difference between the best unrestricted endpoint and reference endpoint is recorded rather than concealed. In both actual summaries the reference is a fixed-coordinate endpoint, with a negligible unrestricted gap:

| Target | Reference unrestricted? | Unrestricted minus reference chi-square | H0 chi-square | H1 baseline chi-square | Full comparison delta chi-square |
|---|---|---:|---:|---:|---:|
| J232128-105122 | False | 1.34666322538e-09 | 75.158887443 | 68.998841878 | 6.160045565 |
| J004131-493611 | False | 1.04852375671e-08 | 95.228460630 | 87.254912074 | 7.973548556 |

These references are numerically equivalent to the recorded unrestricted optima at the stated 1e-4 objective tolerance. They should not be described literally as unrestricted optimizer outputs. A future materially unresolved gap would preclude treating the fixed-coordinate endpoint as an established free optimum. The `profile_baseline_resolved_to_tolerance` flag only checks the objective gap; optimizer success remains a separate requirement.

## Fair H0/H1 restarts

The new H0 starts from the best improved H1 nuisance basin, then H1 starts from the selected H0. This is essential for J004131: the original H0 chi-square is about 110.036, whereas the restarted H0 is about 95.2285. Comparing the old H0 with the new H1 near 87.2549 would exaggerate the improvement. The actual `full_comparison` correctly uses the restarted H0 and gives delta chi-square about 7.97355 for two added shifts. J232128 remains near delta chi-square 6.16005 for four added shifts. These are full-line comparisons, not one-parameter tests of D alone.

A general future-run caveat exists at `common_pair_refine.py:60,104`: `initial_reference` is captured after `refined_h1_from_h0`. If that particular free restart discovers a substantially better new basin and subsequent profile work does not improve it further, the final H0 cross-restart condition will not fire. The actual saved run does not show such a material change: its improved H1 basin was already supplied to the first H0 restart. A generalized implementation should compare against the actual H1 reference used to start the most recent H0, not merely the post-restart reference. No additional optimization is required by this static finding for the current outputs.

## Contour logic and bookkeeping issue

`crossings` computes the first outward crossing on each side of the selected reference, with linear interpolation between actual sampled endpoints. The interpolation formula is correct. It reports brackets and does not claim calibrated coverage. The result describes the connected component around the selected minimum on the sampled grid; it does not enumerate disconnected regions or prove smoothness between points. Levels 1 and 3.84 must remain descriptive objective contours unless coverage is separately calibrated. Active bounds, nonlinearity and the observed change of gas basin invalidate a routine Gaussian/Wilks confidence interpretation.

**Nonblocking bookkeeping ambiguity:** `common_pair_fit.crossings` sorts `bracket_m_s` into ascending coordinate order, but leaves `bracket_delta_chi2` and `fit_files` in inside-to-outside order. On the lower side these orders differ. For example, J232128 level 1 stores coordinates [-497.104360, -474.656673] but delta chi-square [0.753998, 1.023096]; the coordinate-associated deltas are in the reverse order. The computed crossing itself is correct, and file references preserve the association. A consumer must not zip those arrays as corresponding elements. A future schema should either keep all arrays in one common order or expose explicit `inside` and `outside` endpoint records.

The refined J004131 profile uses 21 coordinates, whereas all preserved first and refined fits cover 41 unique coordinates. The older unrefined branches include poor local minima; blindly merging their high objective samples into the newly refined curve creates artificial jumps and spurious first crossings. It is reasonable to preserve them as diagnostics while using the documented refined grid, but the report should say which grid supplies the contours. At common coordinates the source explicitly takes the minimum over all available starts. The final grid has not been optimized at every old off-grid coordinate, so the algorithm does not establish the global profile envelope everywhere.

**Terminology detail:** the `refined_reverse_*` loop traverses coordinates in reverse order but starts each fit from the currently lowest objective endpoint, not from the previous reverse neighboring solution. It is a set of global-best-basin restarts ordered in reverse, rather than a strict reverse continuation or complete hysteresis test. The actual initial parameters are recorded, so this does not affect result reproducibility.

The local unconstrained tangent sigma and constrained nonlinear profile width are different quantities. In particular, the J004131 local sigma near 399 m/s need not match its substantially narrower descriptive profile contour near the active-bound solution. That mismatch should be explained, not averaged away or described as agreement.

## Preservation and resumption

The original summary-export failure is correctly identified: the first-pass `dict(...)` passed `selection` explicitly and again via `**meta`, after the numerical fit products had already been saved. The companion uses a dictionary display with metadata first, removing that duplicate-key call error without editing the historical source or fit products.

The companion writes new `refined_*` products and the refined combined summary. `path_record` changes only loaded dictionaries. The `first_pass_sha256` check protects the first-pass JSON bytes. NPZ preservation is not itself checked by that particular assertion and should be covered by the independent validator. Cached `base.fit` checks source/start hashes, fixed-coordinate identity and NPZ existence. Cached `fit_h0` checks source/start hashes but lacks the analogous NPZ-existence check; all current H0 NPZ products exist, so this is a future interrupted-run robustness detail rather than an observed result defect.

## Source hashes

- `code/common_pair_refine.py`: `f14f4f4b06191c0ca8e2cbaf968862d8202235ff7906cc34cee2906324bebc42`
- `code/common_pair_fit.py`: `d2d8c754dc693d3b61e4acc109fae030fbb3661f9e36269d4329941f8b56b856`

No alpha change, Riemann-clock response, calibrated common-pair significance, or cosmic-time-law fit follows from this source audit.
