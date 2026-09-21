# Internal scientific consistency review of the revised draft

Date: 2026-09-20. Files reviewed: `paper/main.tex`, `paper/theory_sections.tex`, the generated trapped-ion table, `results/analysis_summary.json`, `results/wei2026_summary.json`, `results/numerical_validation.json`, and the associated source/numerical reports. This is a bounded AI-assisted internal review, not external peer review, author approval, or experimental confirmation of the proposed cosmic extension. No manuscript files were edited by this review.

## Conclusion

No critical mathematical, dimensional, numerical, or statistical blocker was found in the reviewed draft. The rewritten argument distinguishes a finite-resource estimation task from a universal physical ceiling, and it does not present the new data as evidence of cosmic aging. Two small wording clarifications were sent to the drafting agent; neither changes the results.

1. The abstract's phrase “covering the first 80 zeros at four drive settings” could imply that all four settings have 80 estimates. The precise distribution is **29 at one setting and 80 at each of three settings**. The main text already gives the correct distribution.
2. The 2026 experiment's dimensionless evolution coordinate is associated with its target Hamiltonian. Its shaped-pulse implementation simulates those target evolutions. The sentence `u = epsilon_0 tau / hbar` should identify the physical-time variable of the target model and should not imply that `u` is itself a measured pulse duration. This is separate from, and consistent with, the manuscript's essential statement that the coordinate is not cosmic age.

## Empirical claims and numerical consistency

The trapped-ion counts are consistent: 29 + 3 × 80 = 269 estimates. The manuscript correctly uses both neighbors to define the final observed index's local gap, including gamma_81 = 202.4935945141405…. This gives a gap of approximately 1.22884257043675 at index 80. The counts of central residuals exceeding a local gap among indices 19–80 are **3, 1, and 2** at drive ratios 8, 12, and 16. These match the current validation output and supersede the earlier audit's truncated-boundary count of two for the first setting.

The generated table agrees with the result CSV. Median absolute errors in the first 18 versus last 20 estimates are approximately 0.1281/0.1280, 0.1398/0.1268, and 0.1914/0.2007. The manuscript properly treats these post hoc groups as descriptive. It neither concludes error stationarity nor attaches a change-point significance to them.

The final-row uncertainties are correctly decoded as 5.76, 0.05, and 0.16. The interpolation uncertainties are not treated as a complete systematic-error budget, and the rounded-zero uncertainty at index 74 is not given infinite weight. The ratio-to-gap line at one is a descriptive reference, not a coverage-calibrated definition of resolved peaks or correct nearest-zero assignment.

The 2026 result consists of 110 processed time-scan points and 18 processed inverse-temperature-scan points. The five reported zero estimates and their residuals agree with the result JSON. Their root uncertainties remain explicitly unavailable. Supplied observable standard deviations are not relabelled as standard errors or root uncertainties. The much higher-index examples remain identified as simulations, not hardware observations.

Index, ordinate, and truncation length are not conflated. The specific values gamma_80, gamma_4200, N(76), and N(4200) agree with the validated outputs. The new reference ordinates are described as high-precision numerical values rather than certified intervals.

## Operational definitions and inference

The pointwise probability criterion specifies a tolerance and an error probability, requires nuisance-parameter coverage, and excludes simply returning a known reference ordinate as an independent analog estimate. Its complete-prefix version separately requires simultaneous coverage under a joint resource budget. Neither individual success nor the end of a scan is silently substituted for a universal ceiling.

The optimized-accessibility monotonicity statement explicitly assumes inclusion of feasible protocol sets together with calibrated performance and fixed targets/criteria. It is correctly presented as conditional set inclusion, not a deduction from cosmic age or thermodynamics alone.

The storage argument distinguishes an explicit table from a compact algorithm. Entropy and operation limits are not converted into a maximum zero index without a task and cost model. The interval-arithmetic RH verification is kept separate from spectroscopy; it supports the narrower rejection of a universal order-4200 limit on physically computable zero information.

## Quantum information and resource envelope

For the stated dimensionless parameter and Hamiltonian `H = hbar kappa theta sigma_z / 2`, `kappa` has inverse-time units and an equatorial pure qubit has QFI `(kappa tau)^2`. With independent trials and the stated local unbiasedness assumptions, the bound `Var(theta_hat) >= 1 / [M (kappa tau)^2]` is correct. The known-visibility quadrature Fisher information is `V^2 (kappa tau)^2` per shot. Unknown calibration, bias, correlated samples, and global aliases are acknowledged rather than suppressed.

The draft explicitly states that the bound is neither attained precision nor a confidence guarantee. The numerical illustration separately **assumes** an attained RMS-like scale `delta^2 = 1/M + b^2`. Its floor amplitudes and normalization are illustrative, not fitted. Although the same letter b is used for a bias in the general MSE model and for a floor amplitude in the illustration, the definitions are local and unambiguous; renaming the illustrative floor would be optional editorial tidying.

Differentiating the smooth zero count gives the stated mean spacing. Solving `delta = q * 2 pi / ln[T/(2 pi)]` gives `T_env = 2 pi exp(2 pi q/delta)`. The draft correctly says this is an asymptotic mean-density crossover, not a local-gap guarantee, bandwidth bound, or attainable experiment. It also states that RMS-to-confidence conversion needs a probability model. Thus the QFI lower bound, illustrative RMS scale, and operational confidence criterion are not equated without notice.

## Published map and cosmic extension

The local published Wang PDF text was checked for the square-root heuristic and spacing claims. The new draft correctly describes the inverse-log-squared schedule as a heuristic ansatz involving an assumed square-root spectral response. It does not promote the numerical correspondence to a strict isomorphism, self-adjoint Hilbert–Pólya construction, or reproduced GUE statistics. The acknowledged poor held-out performance and non-GUE model spacings are consistent with the published source.

The map iteration-to-age identification and the structural-parameter-to-uncertainty response are explicitly additional assumptions. No fitted map coefficient is reused as a physical noise amplitude. The standard-deviation ansatz has a dimensionless logarithm and a separate nonnegative amplitude. Its derivative is correctly a **fractional standard-deviation** rate; the draft does not mistakenly use it as a fractional variance or total-error rate.

The derivative `-p/[t ln(t/t*)]` has inverse-time units. For the illustrative 13.8 Gyr age and Planck-time-scale reference, the quoted rate and half-year change match the independently validated values: approximately -1.0334 × 10^-12 per year and -5.17 × 10^-13, respectively. These are normalized scenarios, not measured drifts.

The MSE decomposition states the zero-mean and mutual-uncorrelatedness assumptions. It acknowledges that an irreducible estimator floor would require a correlation/control model; independent per-shot noise does not remain a floor after unlimited averaging. Single-epoch residuals are not said to identify a cosmic contribution, and publication dates are not treated as acquisition epochs or an aging baseline.

The unrestricted hypothesis permits arbitrarily small amplitude. The draft therefore correctly limits falsifiability to specified positive-amplitude response models; no finite null result is claimed to exclude the entire family. The conclusion does not assert that exponent two is preferred by current data.

## Scope of this review

This check inspected equations and their assumptions, compared the prose with existing independent numeric validation, and consulted the archived 2026 source text and local published Wang text for the specific claims above. It did not independently reproduce the full external experiments, certify all zeta-zero computations, perform an exhaustive literature search, or rerun publication integrity databases. The manuscript retains appropriate author-review draft declarations.

## Drafting disposition

Both wording suggestions were incorporated: the abstract now gives the 29 + 3 × 80 distribution; the NMR paragraph defines tau_model and distinguishes the target evolution coordinate from actual control-pulse duration. No numerical result changed.
