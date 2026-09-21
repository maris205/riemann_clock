# Theory fragment and independent checks

Written 2026-09-20 for `riemann_clock/paper/theory_sections.tex`. This is a new fragment only; original project files were not changed.

## Scope and source provenance

Read the existing `reports/physical_claims_review.md` and `reports/notebook_audit.md`, the published spectral paper's extracted text at `/tmp/riemann_mca_text.txt`, and its bibliographic entry in the previous alpha project. Independently checked the primary metadata/abstracts for [Braunstein and Caves (1994)](https://doi.org/10.1103/PhysRevLett.72.3439), [Giovannetti, Lloyd, and Maccone (2006)](https://arxiv.org/abs/quant-ph/0509179), [Margolus and Levitin (1998)](https://arxiv.org/abs/quant-ph/9710043), and [Platt and Trudgian (2021)](https://arxiv.org/abs/2004.09765). The Ramsey calculation and aging derivative below are elementary derivations for this manuscript; they are not observations from those sources.

Required bibliography keys:

- `Wang2026Spectral`: DOI 10.3390/mca31050193.
- `BraunsteinCaves1994`: DOI 10.1103/PhysRevLett.72.3439, Physical Review Letters 72, 3439–3443.
- `Giovannetti2006`: DOI 10.1103/PhysRevLett.96.010401, Physical Review Letters 96, 010401.
- `Bekenstein1981`: DOI 10.1103/PhysRevD.23.287.
- `MargolusLevitin1998`: DOI 10.1016/S0167-2789(98)00054-2, Physica D 120, 188–195.
- `Lloyd2002`: DOI 10.1103/PhysRevLett.88.237901; author preprint https://arxiv.org/abs/quant-ph/0110141.
- `PlattTrudgian2021`: DOI 10.1112/blms.12460; Bulletin of the London Mathematical Society 53, 792–797.

The fragment uses only standard LaTeX/amsmath/AMS symbol commands and ordinary `\cite`; it introduces no custom macros and has no preamble. All equation and section labels are unique within this fragment. It should be included after the introduction and before, or following, empirical results as the root author chooses.

## Definitions and safeguards

1. `j` is an index, `gamma_j` a dimensionless ordinate, `T` an ordinate threshold, `N(T)` a count, `tau` laboratory interrogation duration, `t` cosmic age, and `n` a map iteration. A finite verified critical-line set is not a proof of RH.
2. A protocol's confidence requirement concerns a signal-coordinate estimator and must include nuisance calibration. A constant estimator returning a published reference value is excluded as independent analog inference. Suggested blinded or held-out validation prevents an exact reference oracle from manufacturing zero measurement error.
3. Isolated-target accessibility and jointly budgeted complete-prefix coverage are different. The prefix requires simultaneous probability control; a pointwise 95% interval at every target is not a 95% guarantee for the entire prefix.
4. Nested **feasible protocols with calibrated performance**, not chronological time by itself, give monotonic optimized accessible height. A larger nominal resource vector is sufficient only if the physical response, control options, and criterion retain the requisite feasibility.
5. Entropy/operation limits do not give a largest index without task complexity. An algorithmic representation can compress a mathematical sequence. The manuscript does not infer a finite mathematical zero sequence.

## Quantum information calculation

For `H = hbar*kappa*theta*sigma_z/2`, an equatorial initial qubit, and interrogation time `tau`, the pure-state QFI is `(kappa*tau)^2`. It is dimensionless per squared dimensionless ordinate; `kappa` has units s⁻¹. Additive QFI for `M` independent trials gives `Var(theta_hat) >= 1/(M*kappa²*tau²)` under regular local-unbiasedness assumptions.

For binary Ramsey outcomes with `p+ = [1+V*cos(phi)]/2`, `phi=kappa*theta*tau+control_phase`, the per-shot Fisher information is

`V²*(kappa*tau)²*sin²(phi) / [1 - V²*cos²(phi)]`.

At quadrature this is `V²*(kappa*tau)²`, hence `sigma >= 1/(V*kappa*tau*sqrt(M))`. Visibility is assumed known and independent of theta. The manuscript does not substitute this lower bound for achieved precision. Global aliases, unknown calibration, biased estimators, and correlated shots are specifically excluded from any universal interpretation.

## Mean-spacing envelope algebra

The derivative of `T/(2pi)*[ln(T/(2pi))-1]+7/8` is `ln(T/(2pi))/(2pi)`. Its reciprocal gives the asymptotic **mean** gap. Solving `delta=q*2pi/ln(T/(2pi))` yields `T_env=2pi*exp(2pi*q/delta)`, with `0<q<1` and delta a dimensionless working resolution scale.

This is only a conditional mean-density crossover. It is not a theorem about local gaps, a fundamental maximum ordinate, a confidence guarantee, or an attained quantum bound. If the main numerical code plots a lower-bound substitution, label it a hypothetical achieved-resolution scenario instead. Bandwidth and signal slopes can restrict observation long before this exponential extrapolation.

## Aging ansatz calculation

The primary ansatz is an **estimator standard-deviation component**:

`delta_c(t)=delta_c0*[ln(t0/t*)/ln(t/t*)]^p`, `t>t*`, `delta_c0>=0`, primary `p=2`.

For positive amplitude, `d ln(delta_c)/dt=-p/[t ln(t/t*)]`. The derivative of the corresponding **variance** is twice this fractional rate. This distinction is essential when fitting a variance rather than a standard deviation.

Using `year = 365.25*86400 s`, `t0=13.8e9 year`, and illustrative `t*=5.391247e-44 s` gives:

| Quantity | Independently calculated value |
|---|---:|
| `t0` in seconds | `4.3549488e17` |
| `ln(t0/t*)` | `140.24422681374799` |
| Present fractional SD derivative, p=2, per year | `-1.0333939551348218e-12` |
| Half-year linearized relative SD change | `-5.166969775674109e-13` |
| Stable exact log SD ratio over half a year | `-5.166969775579837e-13` |
| `delta_c(1 Gyr)/delta_c(t0)` | `1.0385075690587888` |
| `delta_c(100 Gyr)/delta_c(t0)` | `0.9723435970864454` |

The stable finite-interval calculation uses

`-p*log1p(log1p(dt/t0)/ln(t0/t*))`

for `ln[delta_c(t0+dt)/delta_c(t0)]`, avoiding cancellation at tiny `dt/t0`. For the relative change itself use `expm1` of that result. These are illustrative values, not a Planck-scale derivation or measured drift.

The MSE decomposition is stated only for zero-mean mutually uncorrelated random components plus a remaining bias squared. A floor needs temporal/cross-shot correlations or another justified mechanism that prevents averaging it down. Independent per-shot cosmic noise is not an irreducible estimator floor. No noise correlation spectrum is invented in the manuscript.

## Connection to the published map and identifiability

The published map schedule is a heuristic numerical construction, with a square-root spectral response assumed in its motivation. The paper explicitly disclaims mathematical isomorphism and reports local spacings that do not reproduce GUE. The new fragment therefore claims only **motivation for testing a functional form**. Mapping map index to physical time and map parameter to a standard deviation are two additional assumptions; even retaining the same exponent requires a specified response relation. Map-fit coefficients are not copied into a dimensioned experimental floor.

One-epoch residuals cannot isolate a common cosmic contribution from unknown calibration/model discrepancy. They cannot choose exponent two. Different papers appearing six months apart are not controlled observations of one device at different cosmic ages. A positive-amplitude, fully specified model can be excluded; the family permitting amplitude arbitrarily close to zero is not falsified by any finite null result.

No new threshold near 4200, next-zero waiting time, or universal cutoff is predicted. All comments about the old 4200 value are diagnostic exclusions rather than alternative forecasts.
