# Audit of the two original Riemann-zero notebooks

Audit date: 2026-09-20. The original notebooks and PDF were read without modification. Neither notebook was executed. Small independent calculations below use only the literal table, ordinary arithmetic, and three `mpmath.zetazero` reference values. This report audits the old implementation; external verification of the hardcoded table must be supplied separately.

## Files and provenance

- `Cosmic-Chaos-Alpha/4-riemann_error_evolution_nyquist.ipynb`: SHA-256 `03a7ded20152c31ffe4991916bad05fd1269c2b5fb50e71c656b3d4557eb972d`.
- `Cosmic-Chaos-Alpha/5-theoretical_observable_height_heatmap.ipynb`: SHA-256 `b2e3e6e9f1af7b23450e70ba014135ef7d8e21e6a4681db2d78da4a88fab565b`.
- `Cosmic-Chaos-Alpha/theory_result/generate_figures.py`: SHA-256 `432929824dd3fdd23b26def94056e83426172fb47d3563e3aaf3dcbcc7118b4a`.

The actual two notebook names begin with `4-` and `5-`. Each has one populated code cell and an empty second cell.

## Error notebook: what the code actually does

The first notebook contains a literal CSV-style table with 80 zero indices and reference ordinates rounded to approximately three decimals. It includes four experimental-looking columns, of which the first is populated only for indices 1–29. The plotted columns are labelled drive settings Ω=8, 12 and 16. The notebook does not identify the source table, DOI, acquisition dates, units for Ω, or meaning of the parenthetical uncertainties. Their experimental provenance cannot be established merely from the hardcoded numbers. This literal table is not generated randomly within this notebook.

Parenthetical uncertainties are parsed in the conventional last-digit notation. In particular `199.24(576)` means a central value of 199.24 and an uncertainty of 5.76, **not** an uncertainty of 576. Conversely, Ω=8 at index 74 is `189.50(0)` and becomes an exact zero uncertainty in the parser; this must not become an infinite inverse-variance weight. The original source must establish whether this is a rounding artifact or transcription error.

The code computes central residuals (experimental central value minus rounded mathematical ordinate) and draws uncertainty bars. It performs no parameter fitting, covariance estimation, break-point test, likelihood comparison, uncertainty calibration, or significance calculation. The cutoff index 18.5 is entered manually. Everything to its right is manually shaded and called “Phase Aliasing & Breakdown Region (Error > Zero Spacing)”. The comments and annotation calling index 19 “Early Topological Breaking”, and index 80 “Complete Decoherence”, are interpretations inserted in plotting code, not inferences produced by an analysis.

The fixed y-axis range −4.5 to +4.5 also clips large uncertainty bars. A revision should expose rather than clip unusually uncertain measurements, for example in a separate uncertainty-to-spacing panel.

### Transparent descriptive recomputation

For each point, define the comparison spacing as the smaller of the adjacent spacings in the rounded reference table. At the two boundaries, use the single available adjacent spacing. This is a descriptive neighboring-gap comparison, not a detection or resolvability theorem, and does not use uncertain estimates as though they were independent Gaussian observations.

| Setting | Median absolute residual, indices 1–18 | Median absolute residual, indices 61–80 | Median reported uncertainty, indices 1–18 | Median reported uncertainty, indices 61–80 |
|---|---:|---:|---:|---:|
| Ω=8 | 0.1285 | 0.1280 | 0.030 | 0.195 |
| Ω=12 | 0.1400 | 0.1265 | 0.060 | 0.180 |
| Ω=16 | 0.1915 | 0.2005 | 0.065 | 0.375 |

The median central residual does not show a universal high-index cliff. Some reported uncertainties become larger; that observation alone does not identify an underlying cosmic effect or a hard cutoff.

Of the 62 measurements beyond index 18, the absolute central residual exceeds the smaller adjacent reference spacing for only 2 points at Ω=8, 1 point at Ω=12, and 2 points at Ω=16. The corresponding reported uncertainties exceed that spacing for 10, 12, and 12 points. Thus the blanket shaded region labelled “Error > Zero Spacing” is not justified under either reading.

Across all 80 measurements, the median absolute residuals are 0.1405, 0.1695, and 0.1955; root-mean-square residuals are 0.41552, 0.36805, and 0.35625 for Ω=8, 12, and 16, respectively. These three modes are observations of the same zero sequence under related experimental controls. They must not automatically be treated as independent cosmological experiments.

## Heatmap notebook: an assumed expression, not a numerical derivation

The only mathematical calculation is

\[
T(X,\epsilon)=\frac{\ln X}{\pi\sqrt{\epsilon}}.
\]

The notebook labels this a purely physical derivation, but it provides no sampling protocol, derivation, noise model, explicit-formula evaluation, prime computation, comparison likelihood, or cosmic-time calculation. It evaluates the expression on a 200×200 grid, with `X` from 10^10 to 10^80 and ε from 10^−5 to 10^−1, then clips values above 15,000. Approximately 0.585% of the grid values are clipped.

The two labelled anchors are chosen inputs:

- `X=10^10`, ε=10^−2 gives `T=73.2935599`.
- `X=10^60`, ε=10^−4 gives `T=4397.6135933`.

The laboratory marker is actually drawn at `X=10^11`, despite its label saying `10^10`; the plotted coordinate corresponds to `T=80.6229159`. This displacement is cosmetic and should not be reproduced as a data point.

No apparatus quantity is measured and identified with the prime upper bound X in the notebook. Likewise, no noise calibration produces ε. Neither the scale `10^60` nor an α-variation amplitude is a demonstrated universal physical precision of the universe.

Nyquist sampling requires a specified independent variable, sample spacing, band limit, and reconstruction class. For a prime-related Fourier representation with phase `T log p`, the logarithmic window extent and the sampling mesh play different roles. The window length governs Fourier resolution, while the mesh governs aliasing. The displayed formula does not follow from invoking the name “Nyquist” without those additional definitions and assumptions.

## Inconsistencies with the old PDF

1. **Index and height are different.** Independently computed ordinates are `γ18=72.0671576745`, `γ19=75.7046906991`, and `γ80=201.2647519437`. Therefore a cutoff at ordinate 73.3 lies between indices 18 and 19, not near index 80. The error notebook maps these correctly in its vertical-line position; the PDF compares a purported ordinate near 76 with the index 80 and treats their numerical proximity as corroboration.
2. **The numbers 4200 and 4397 come from different prescriptions.** The PDF uses `ε ∝ T^(−a)` with a chosen `a≈1.15` and benchmark `T=4200`, ε=10^−4. Reducing precision to ε=10^−2 yields `T=76.5806160`. An exponent of `1.1626830` would map the selected pair exactly to 80. There is no exponent estimation or uncertainty calculation in either notebook. The heatmap instead assumes the logarithmic-X expression above, producing 4397.6 rather than 4200.
3. **A count near 4200 is not an ordinate near 4200.** The smooth Riemann–von Mangoldt count `T/(2π)[ln(T/(2π))−1]+7/8` is approximately 3680.67 at `T=4200`. Its value is an asymptotic smooth approximation, not the exact integer count.
4. **The claimed 100-million-year wait is not reproduced.** Under the PDF's explicit assumptions `X∝T²` and `X∝cosmic age`, a current age of 13.8 Gyr gives `13.8 Gyr × [(4201/4200)²−1] = 6.5722 Myr`, not 100 Myr. No notebook performs that time calculation. A change in ordinate by one must also not be described automatically as access to one additional zero.
5. **The far-future numerical claim is incompatible with the heatmap law.** At fixed ε=10^−4, changing X from 10^60 to 10^120 doubles the heatmap height to 8795.23; it does not give 10^5. The PDF and code need one declared conditional scenario, not mutually different scaling laws.
6. **The explicit-formula argument changes the meaning of its bound.** A truncation depth sufficient to approximate a prime-counting function at a selected X and tolerance is a computational accuracy requirement. It is not by itself a maximum physical height of observable zeros. The assertion of exponential convergence and a universal unique cutoff is not implemented or established by the code.
7. **The main experimental citation was unfinished.** The PDF's reference [6] is explicitly a placeholder for a trapped-ion paper. Before any new inference, the actual source and supplementary table must be verified.

## Separate generated figures contain simulated, not measured, errors

`theory_result/generate_figures.py` explicitly constructs “USTC ion-trap error data” using pseudorandom draws with seeds 10, 20 and 30. The code prescribes three regions: small noise for n<19, increasing scatter for 19≤n<50, and “cliff-like divergence” for n≥50. Its uncertainty bars are also prescribed piecewise functions. The `omega` argument is unused inside the generator; the three sets differ through their seeds, rather than a physical drive-frequency response.

Those synthetic arrays are then presented in figures titled “USTC Ion-Trap Error Evolution” or “USTC Ion Trap vs. Nyquist Phase Truncation”. These figures cannot be used as experimental verification. The simulated cliff is a construction of the generator. Do not conflate these outputs with the literal table in notebook 4. A new revision should rebuild every experimental figure directly from a provenance-verified table; any retained simulations must be clearly labelled synthetic demonstrations.

## Defensible revision direction

Preserve the exploratory question: under a declared physical measurement protocol, how does a finite resource/noise budget limit the precision or range of estimates of mathematical zeros, and how would an assumed cosmic-age dependence of those resources change that range? Keep the exact mathematical zero sequence independent of the observer. Distinguish experimental control limitations, statistical uncertainty, systematic bias, and speculative cosmic evolution.

Use the new published dynamical-systems paper as motivation for a conditional nonautonomous noise or structural-parameter model. It does not supply a proof of a universal cosmic cutoff, an apparatus-independent Nyquist limit, or a direct identification of its iteration index with cosmic time. Report uncertainty growth in the ion experiment descriptively unless the original measurement model supports a stronger inference. The present notebooks do not test temporal evolution, because their rows vary in zero index and drive setting rather than cosmic observation epoch.

## Subsequent source verification and endpoint convention

The primary supplement was subsequently verified and transcribed; see `experimental_sources_audit.md`. The final analysis uses independent high-precision reference ordinates, including gamma_81 to define the forward gap at index 80. Its beyond-18 central-error counts are consequently **3, 1, 2**, whereas the original notebook-only comparison above used a one-sided terminal gap and gave **2, 1, 2**. The final manuscript and `numerical_validation.md` use the former convention. Original rounded reference-column discrepancies at indices 22 and 23 are retained in the source CSV and documented.
