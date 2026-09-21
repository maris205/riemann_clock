# Independent numerical and interpretation check: resource heatmap

Checked on 2026-09-20. This check concerns an illustrative resource model, not a fit to either experimental dataset.

## Formula and meaning

The existing manuscript defines the asymptotic mean spacing as `dbar(T) = 2*pi / ln(T/(2*pi))`. Equating a stipulated achieved dimensionless resolution `delta` to `q*dbar(T)` gives

```text
delta(M,b) = sqrt(1/M + b^2)
log10(T_env) = log10(2*pi) + 2*pi*q / (delta*ln(10))
q = 0.25
```

This algebra is correct. The color encodes the logarithm of an ordinate envelope, not a zero index, a zero count, or an experimentally established maximum. The resolution model is assumed attained; its use does not turn a Cramer–Rao lower bound into an attainment guarantee.

Increasing `M` at fixed positive `b` decreases `delta` and increases the envelope toward a finite plateau. Lowering `b` at fixed `M` increases the envelope. Equal statistical and floor contributions occur at `M = 1/b^2`. The plateau is

```text
lim[M -> infinity] log10(T_env) = log10(2*pi) + 2*pi*q/(b*ln(10)).
```

These statements assume the stipulated floor survives the repetitions under consideration. A per-shot independent random error would not in general satisfy that assumption.

## Independent arithmetic

Values below were evaluated using Python's standard-library `math` functions, independently of the plotting implementation.

| Location | M | b | delta | log10(T_env) |
| --- | ---: | ---: | ---: | ---: |
| Minimum over displayed domain | 1 | 1 | 1.414213562373095 | 1.280559754304186 |
| Maximum over displayed domain | 10000 | 0.03 | 0.031622776601684 | 22.370864187438325 |
| A | 100 | 0.3 | 0.316227766016838 | 2.955448300266136 |
| B | 3000 | 0.3 | 0.300555042102663 | 3.067941081307259 |
| C | 3000 | 0.05 | 0.053229064742238 | 13.614264092727336 |

The full domain is `1 <= M <= 10^4`, `0.03 <= b <= 1`. Monotonicity places the exact extrema at the stated corners; the color limits should retain both extrema without clipping. The plateau at `b=0.3` is `log10(T_env)=3.072140458094517`, explaining why increasing repetitions from A to B yields modest improvement. The plateau at `b=0.05` is `14.441943406776526`.

## Interpretation requirements

- A, B, and C are invented resource settings for explanation and must not be described as experimental or cosmological anchors.
- The mean-spacing formula does not impose local gaps, signal strength, bandwidth, computational preparation cost, calibration, or global identifiability. The large upper-envelope values are an illustration of its exponential sensitivity, not a prediction of reachable zero height.
- At the lowest corner, the formal envelope is only `T_env ~= 19.08`; the asymptotic spacing approximation is not quantitatively strong there. The map should be labeled an algebraic illustration using asymptotic mean density throughout.
- There is no cosmic-age coordinate or cosmic-aging fit in this heatmap. It neither establishes a cosmic measurement floor nor tests the proposed `1/ln^2(t/t_*)` schedule.
- Negative or zero floors are excluded from this particular logarithmic-axis map. No claim about a universal positive floor follows from that plotting choice.

## Implementation inspection

Inspected `code/plot_resource_heatmap.py`, `results/resource_heatmap_grid.npz`, and `results/resource_heatmap_summary.json` after generation. The implementation passes the numerical check:

- Independently recomputed all 77,361 grid entries using scalar standard-library arithmetic. The maximum absolute difference in `log10(T_env)` was `7.11e-15`, consistent with floating-point roundoff.
- Verified shape `(241, 321)`, the intended positive domains, strict increase with repetitions, and strict decrease with assumed floor.
- Verified both full-domain extrema match the summary exactly. The plot uses those extrema as `vmin` and `vmax`; no grid truncation or color-range clipping is introduced.
- Verified all three star values independently and confirmed their coordinates are strictly inside the plotted domain.
- The axis specifies ordinate units, the color label specifies `log10(T_env)`, the subtitle identifies the conditional model, and the stars are explicitly illustrative. The output summary correctly identifies the quantity as neither an observed cutoff nor a zero count.
- The dashed line `b=M^(-1/2)` correctly marks equal contributions to the assumed error budget.

The only optional visual suggestion was to move the start of the text "Precision floor dominates" slightly rightward so that its left edge also lies on the floor-dominated side of the dashed crossover. Root handles final visual and PDF inspection. The low-height asymptotic limitation described above remains a matter for the caption/text, not a numerical implementation failure.
