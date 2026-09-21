# J064326−504112: amended failure-case displacement diagnostic

**Adding three relative line shifts gives a limited objective improvement and does not repair the weak 1611 mismatch. The target remains outside the precision/common-displacement sample.**

This is a root-authorized protocol amendment after the initial quality screen. The earlier adequacy gate decided whether a model was suitable for a precision displacement estimate; it was not a population-level null test. Omitting every inadequate H0 without checking H1 could hide genuine sensitivity to frequency changes. The present bounded comparison explicitly tests that possibility for this target and retains the result, without promoting it to a calibrated measurement.

The amendment was recorded before the new fits in `protocol_amendment.json`. It preserves the original four lines (1608, 1611, 2374, 2382), all 436 native pixels, the [−145,+100] km/s window, 16 shared gas components, original expected-fluctuation diagonal weights and 21 subpixel integration. The weak 1611 line is retained. There is no residual clipping or component reselection.

| Model | Parameters | Pixels | χ² | Nominal dof | Weak 1611 χ²/pixel |
|---|---:|---:|---:|---:|---:|
| Selected original H0 | 64 | 436 | 532.225273754 | 372 | 2.058481575 |
| H1, three relative shifts | 67 | 436 | 527.587722694 | 369 | 2.031517966 |

The improvement is Δχ² = **4.637551060 for 3 added parameters**. Both models still fail the same weak-line residual gate (χ²/pixel≤1.8). Allowing frequency offsets therefore does not repair this model failure within the frozen architecture and bounds. This is not an exhaustive test of every possible frequency-dependent effect, blended-line model or gas architecture.

The diagnostic H1 offsets relative to fixed 2374 are:

| Transition | Conditional best-fit offset (m/s) |
|---|---:|
| 1608 | +5.548692 |
| 1611 | -895.272585 |
| 2382 | -217.020141 |

These are nuisance-fit endpoints under an inadequate model. They have no reported precision uncertainty and are not alpha estimates or cosmic-age observations. In particular, the 1611 shift is close to the −1000 m/s permitted edge, though it is not classified as an active shift bound. The largest localized weak-line residual remains conceptually different from shifting the entire weak profile.

Only two new fits were run: H1 from the original primary H0 (45 evaluations, ftol stop), followed by H0 cross-started from H1 (36 evaluations, ftol stop). The cross-started H0 has χ²=532.225296838 and does not beat the original H0; the predeclared optional extra H1 was therefore not triggered.

H0 and H1 both retain active `logb_14` (lower 0.5 km/s) and `zero_1611` (upper +0.02). Their recorded first-order optimalities are 0.210247 and 0.144307. Successful ftol stopping does not establish strong stationarity or a global optimum. 16 components is the upper end of the previously searched count grid.

All original 20 fit JSON/NPZ files and the frozen primary configuration remain protected. The original gated-run reports describe that historical protocol; this separate amendment supplies the additional explanatory H1 result. The canonical primary status remains `not_estimated_primary_due_model_failure`. The present outcome cannot be described as a population null result.

Reproducibility: `summary.json` points to the selected original H0 and new H1; new fit records and arrays are under `J064326-504112/`. `failure_pair_comparison.pdf` shows the same raw pixels, both models and residuals.

Independent review passes 93/93 grouped checks, including preservation of all 40 original fit files, objective and Jacobian replay, and 28 representative finite-difference columns including all three shifts. See `independent_review.md` and `independent_review.json`. Calculation verification does not establish physical validity.
