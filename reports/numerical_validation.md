# Independent internal numerical validation

Date: 2026-09-20. Scope: the descriptive numerical products in this project, the downloaded public experimental summaries, and the conditional illustration formulas. This is an independent internal consistency check; it is not external peer review, an author approval, or a physical validation of the cosmic-time hypothesis.

The numerical reanalysis passed the checks below. The analysis does not fit a cosmological parameter, apply a joint Gaussian likelihood to the quoted interpolation uncertainties, or infer a universal physical cutoff.

## Reproducible check

Run from this directory's parent:

```sh
python code/validate_analysis.py
```

The validator reads existing products without importing either analysis program. It independently parses the supplementary PDF with decimal arithmetic, recalculates all mathematical references at 50 decimal digits, checks the NPZ data directly with `allow_pickle=False`, and writes `results/numerical_validation.json`. The source programs remain unchanged by this review.

## Input integrity and transcription

- All three original PDF/notebook inputs match both SHA-256 and byte counts in `data/source_project/source_manifest.json`. The original files under `Cosmic-Chaos-Alpha` remain intact.
- All 14 entries in `data/raw/experimental_sources_download_manifest.json` match their recorded hashes and byte counts. This statement concerns the listed entries, not every file present in the directory.
- Direct parsing of Tables I and II of the He et al. supplementary PDF reproduces all **269** extracted records: 29 at drive ratio 5 and 80 each at ratios 8, 12, and 16.
- The parenthetical uncertainties in the last row are **5.76, 0.05, and 0.16**, respectively. In particular, `199.24(576)` denotes an uncertainty of 5.76 in the reported ordinate, not 576.
- The sole uncertainty rounded to zero is `189.50(0)` at index 74 and drive ratio 8. It is flagged and receives no infinite inverse-variance weight. Its uncertainty marker is explicitly omitted from the logarithmic ratio panel.
- Wei et al.'s two processed time scans contain 110 unique rows. All **880** numeric cells (beta, encoded time, and six means/standard deviations) exactly match their NPZ inputs. The processed beta scan contains 18 rows and all **144** numeric cells exactly match its NPZ input. No upstream Python source was executed.
- The five reported 2026 root estimates are 14.12, 20.96, 25.09, 30.44, and 32.93; an additional independent source-text check confirmed that these agree with the archived publisher Results paragraph. Their root uncertainties are unreported and remain explicitly missing. Measurement error bars from the scan arrays are not silently converted into fitted-root uncertainties.

## Mathematical references and source-table corrections

All 83 stored ordinates (indices 1–81, 4200, and 4201) were recalculated at 50 decimal digits. The largest absolute difference from the stored 35-digit values was approximately **3.80 × 10⁻³²**. The largest evaluated residual `|zeta(1/2 + i gamma_j)|` was approximately **1.50 × 10⁻⁴⁷**. These are numerical checks using the same underlying `mpmath` library at higher precision; they are not certified interval enclosures or independent software proofs of completeness.

The original supplement's rounded reference column has two entries outside half a unit in the last printed decimal place:

| Zero index | Printed reference | Recomputed ordinate |
|---|---:|---:|
| 22 | 82.914 | 82.9103808540860… |
| 23 | 84.736 | 84.7354929805171… |

Index 23 is only about 0.000007 beyond the nominal rounding tolerance; no experimental significance should be attached to this small discrepancy. The source transcription preserves both printed entries. New residuals and neighboring gaps use the independently recomputed ordinates, which is the appropriate separation between provenance and analysis.

The count/ordinate distinction is numerically verified:

- `gamma_80 = 201.2647519437038…`, not 80.
- `gamma_4200 = 4697.3541474715…`, not 4200.
- The reference values bracket 76 between ordinates 19 and 20, and 4200 between ordinates 3681 and 3682, agreeing with `N(76) = 19` and `N(4200) = 3681`.

## Residual and uncertainty diagnostics

All 269 derived residuals, local gaps, and normalized ratios match an independent calculation. All 15 group-summary rows match independent medians and direct counts.

| Drive ratio | Median absolute residual, indices 1–18 | Median absolute residual, indices 61–80 |
|---|---:|---:|
| 8 | 0.128148 | 0.127996 |
| 12 | 0.139788 | 0.126783 |
| 16 | 0.191433 | 0.200734 |

These descriptive values do not display the claimed universal breakdown immediately after index 18. They do not prove that errors are stationary or that no apparatus-dependent degradation exists. The reported interpolation errors and the central residuals are different diagnostics; one cannot substitute for the other.

The local gap is the smaller of the distances to both adjacent mathematical zeros, with a one-sided definition only at index 1. Crucially, the final observed index 80 still has a mathematical neighbor at index 81:

```text
gamma_81 = 202.4935945141405…
min(gamma_80 - gamma_79, gamma_81 - gamma_80) = 1.22884257043675…
```

For indices 19–80, central residuals exceeding one such gap occur at:

- Drive 8: indices **41, 79, 80**, or **3 of 62**.
- Drive 12: index **54**, or **1 of 62**.
- Drive 16: indices **35, 72**, or **2 of 62**.

The first notebook audit reported 2 instead of 3 for drive 8 because it used only the backward gap at the last index of its truncated reference table. Using the next mathematical zero corrects that boundary convention. The new manuscript and numerical products should use **3, 1, 2**. This small change does not alter the main finding that the data do not establish a universal cutoff.

The plotting code permits the complete vertical extent of the original error bars. It does not impose the old “breakdown” shading or the old clipped ordinate range. The ratio panels use logarithmic axes and omit only the explicitly flagged zero uncertainty marker.

## Conditional resource illustration

All 483 scenario rows obey the declared precision assumption

```text
delta² = 1/M + b²
```

and the inverse mean-spacing relation

```text
delta = q [2 pi / ln(T_env / (2 pi))],  q = 0.25.
```

The plotted expression for `log10(T_env)` is algebraically correct. Its values increase with independent repetitions and approach a finite limiting value when the assumed floor is positive. The choice `kappa * tau * visibility = 1`, the statistical attainment assumption, floor values, and tolerance q are scenario inputs. The lower Cramér–Rao bound alone does not establish that any real protocol attains these values. The relation describes the mean density of zeros, not a guaranteed resolution of every local gap, a universal maximum ordinate, or an experimentally calibrated cosmic bound.

## Conditional cosmic-time illustration

All 603 normalized aging rows agree with direct high-precision evaluation for exponents 1, 2, and 3. With the specified illustrative choices `t0 = 13.8 Gyr`, `tstar = 5.391247e-44 seconds`, and exponent 2, independent numerical differentiation gives

```text
d ln(delta_cos) / dt = -1.03339395513482e-12 per year.
```

Direct 50-digit evaluation at `t0 + 0.5 years`, avoiding double-precision cancellation, gives

```text
delta_cos(t0 + 0.5 yr) / delta_cos(t0) - 1
  = -5.16696977557850e-13.
```

The signs, time-unit conversion, and sub-year evaluation are correct. The normalized illustration does not determine a floor amplitude, identify the physical clock with a dynamical-map iteration, or measure an actual cosmic aging effect. The data do not select exponent 2 over its alternatives.

## Review conclusion

No numerical or transcription bug was found in the checked current analysis products. The consequential corrections are scientific interpretation and provenance: maintain the distinction between zero index and ordinate; retain the source-table discrepancy note; report the final-index gap convention correctly; distinguish measurement uncertainty from fitted-root uncertainty; and label all resource/aging curves as conditional scenarios. The preserved old waiting-time numbers are arithmetic consistency diagnostics under discarded assumptions, not new physical forecasts.
