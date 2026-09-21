# 2026 NMR experimental figure and root residual table

All outputs in this note are confined to `riemann_clock`. The original project and the downloaded upstream scripts were not modified. No upstream script was executed.

## Source and scope

The source is Shijie Wei et al., “The Riemann Hypothesis manifested in dynamical quantum phase transitions,” *Nature Communications* **17**, 8163 (2026), [DOI 10.1038/s41467-026-74935-8](https://www.nature.com/articles/s41467-026-74935-8). Formal publication was **1 July 2026**, within the six-month window ending 20 September 2026. Its [arXiv preprint 2511.11199](https://arxiv.org/abs/2511.11199) first appeared **14 November 2025**, outside that window. This is a recent formal publication, not a claim that the results first became public in the last six months.

The physical experiment uses five NMR qubits and examines the first five zero neighborhoods. The discussion of the trillionth zero in the source concerns numerical simulation. The two are not interchangeable experimental records or evidence for an increase in the largest physically observed index.

The authors' [public data repository](https://github.com/lqf2025/Riemann-data) is pinned to commit `a8b0b38202d6f330297212f8755c8df04242f9d6`. The retained input products are explicitly classified by its README as **processed experimental data**. They contain means and standard deviations, not the original repeated NMR signals or a full covariance matrix.

## Figure construction

Run from any working directory:

```bash
python /root/autodl-tmp/dsc_world/riemann_clock/code/plot_recent_experiment.py
```

Dependencies are Python, NumPy, mpmath and Matplotlib. The script uses the noninteractive Agg backend and makes no network requests. `code/extract_public_experiments.py` reconstructs its CSV inputs from the retained publisher supplement and author NPZ products.

Outputs are `figures/wei2026_experiment.pdf`, `figures/wei2026_experiment.png`, `results/wei2026_summary.json`, `results/wei2026_published_root_residuals.csv` and `results/wei2026_root_table.tex`. The TeX file contains a complete `tabular` environment using `booktabs`, so it can be included inside a table float without relying on fragile row-only input.

Panels (a)–(e) show all **110** processed time-scan points: five neighborhoods, 11 points per neighborhood, at each of beta = 0.3 and 0.5. The ordinate is the authors' `abs_mean`, labeled probe coherence magnitude, with no rescaling or reconstruction from separately averaged x/y components. The bars are `abs_std`; the source's Figure 4 caption identifies these as **one standard deviation from repeated experiments**. They are not standard errors, confidence intervals for roots, or a complete systematic error budget. No smoothing, connecting fit, interpolation or root extraction was performed.

The horizontal coordinate is the authors' encoded dimensionless time minus an independently computed reference ordinate, u − γⱼ, with u = ε₀τ/ℏ and laboratory duration τ. It is **not cosmic age**. The reference values come from `mpmath.zetazero` at 45 decimal digits. Grid centers differ from these values by less than 10⁻⁶, as expected from the supplied six-decimal centers. Dashed vertical lines mark known mathematical reference values rather than fitted measurements. The points themselves retain their original horizontal positions; they are not artificially shifted by beta to separate the markers.

Panel (f) shows all **18** processed beta-scan points, at fixed u = 14.134725142. Its vertical scale differs from the first five panels and is fully labeled. This is a separate scan. In particular, its value at beta = 0.5 should not be pooled with the corresponding central time-scan value as though the two were identically conditioned repetitions. The downloaded subset does not supply the covariance or run-level calibration needed for such an analysis. Both are plotted in their reported normalization.

The first-five coherence minima are nonzero at the supplied sample points. Finite encoding, sampling, calibration and noise prevent an inference that these sample values are estimates of a universal physical floor. Distinguishing off-critical-line beta from critical-line beta across the whole designed response is also a different question from testing a cosmic-age-dependent resolution law. The figure illustrates measured protocol responses; it does not fit cosmic evolution.

## Published root estimates against independent references

The source reports polynomial-fit root estimates rounded to two decimals. The script **transcribes** those published values, subtracts high-precision mathematical references, and supplies no invented root standard deviations. The calculated residuals are descriptive and are not a new root fit.

| Index j | Mathematical reference gamma_j | Published estimate | Estimate minus reference |
|---:|---:|---:|---:|
| 1 | 14.134725142 | 14.12 | −0.014725142 |
| 2 | 21.022039639 | 20.96 | −0.062039639 |
| 3 | 25.010857580 | 25.09 | +0.079142420 |
| 4 | 30.424876126 | 30.44 | +0.015123874 |
| 5 | 32.935061588 | 32.93 | −0.005061588 |

The mean absolute difference is approximately **0.035219**, and the largest absolute difference is approximately **0.079142**. The estimates' two-decimal publication rounding is retained; extra digits in residuals serve reproducibility and do not imply measurement accuracy. Root uncertainties are not provided for these five estimates, so neither a chi-squared fit nor a significance claim is justified by this table alone. Different encodings, protocols and uncertainty conventions also prevent a direct accuracy ranking against the 2021 trapped-ion measurements.

## Validation performed

The script checks finite data, nonnegative published standard deviations, the exact 110/18/5 input counts, both beta groups, all 11 unique point indices in each neighborhood, time increments of 0.3, and agreement of the grid centers with independent reference zeros. It preserves SHA256 values of the three input CSVs, the source commit and the mpmath version in the JSON output. Both PDF and PNG were generated successfully and the PNG was visually inspected for legibility, full error bars and panel labels. No point was excluded.

Suggested manuscript caption:

> Processed NMR measurements from Wei et al. (2026). Panels (a)–(e) show the authors' reported probe-coherence magnitudes in the first five zero neighborhoods for beta = 0.5 (circles) and beta = 0.3 (squares), plotted against encoded-time offset from the mathematical ordinate. Panel (f), with a different vertical scale, shows the separately acquired beta scan near the first zero. Bars are the published one-standard-deviation values from repeated experiments. Dashed lines identify mathematical reference parameters. These are 128 processed measurement points in total; no additional root fitting was performed. The dimensionless encoded time is distinct from cosmic age, and the source's high-index simulations are not included as experimental measurements.
