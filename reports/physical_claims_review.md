# Physical interpretation and recent primary literature review

Audit date: 2026-09-20. Scope: the original `riemann_zeros_Fixed.pdf`, using the preflight-passed text extracted to `/tmp/riemann_zeros_original.txt`; independent review of the physical claims and targeted primary-source literature search. This review does not establish an exhaustive systematic review or a new cosmological measurement.

## Recommended central claim

**A fixed physical protocol, with finite duration, signal strength, calibration accuracy, and sampling resources, estimates encoded Riemann-zero ordinates with finite uncertainty. A protocol-dependent accessible height can increase when these resources improve. A link between that resource growth and cosmic age is an additional, presently uncalibrated hypothesis.**

This preserves the motivating idea while separating three assertions:

1. Mathematical zeros and their ordinate values are time-independent objects.
2. A laboratory estimator of a signal encoding those values has an uncertainty budget.
3. A proposed age dependence of the physical resource budget is a testable extension, not a result of items 1 and 2.

Finite age alone does not supply an established nonzero intrinsic error floor for every conceivable protocol. Nor does a failed estimate make a mathematical zero nonexistent. The phrase “not accurately observable” must mean failure of a declared tolerance and confidence criterion, rather than lack of an exact infinite decimal expansion: no ordinary measurement reports such an expansion in the first place.

## Essential notation and internal numerical checks

Write a zero as `rho_j = beta_j + i gamma_j`. `j` is an index; `T` is an ordinate threshold; `N(T)` counts zeros with positive ordinates at most `T`, with multiplicity. Riemann–von Mangoldt gives

\[
N(T)=\frac{T}{2\pi}\log\frac{T}{2\pi}-\frac{T}{2\pi}+O(\log T).
\]

The smooth count has the additional constant `7/8`; neither formula licenses `T ≈ j` for the original comparisons. At a zero, endpoint conventions must be stated when using a smooth formula plus the argument term.

Independent numerical checks using `mpmath` at 30 decimal digits:

| Quantity | Value |
|---|---:|
| `gamma_80` | 201.264751943703788733016133428 |
| `gamma_4200` | 4697.3541474715424375704997813 |
| `gamma_4201` | 4697.77440376115430302225894068 |
| `N(4200)` | 3681 |
| `N(76)` | 19 |

These are high-precision numerical checks, not newly certified interval enclosures. In particular an ordinate threshold of 76 cannot be compared directly with observation of the first 80 zeros.

The old waiting-time claim also fails its own elementary scaling. If `X ∝ age` and `X ∝ T²`, taking `age_0 = 13.8 Gyr` and changing `T` from 4200 to 4201 yields

\[
\Delta t=t_0[(4201/4200)^2-1]=6.5722109\ \mathrm{Myr},
\]

not 100 Myr. Using the actual 4200th and 4201st ordinates gives 2.4693885 Myr. These calculations expose inconsistent substitutions; **neither is a physical forecast**, since the required time and resource mapping has not been derived or measured.

## An operational definition that can actually be tested

For each protocol `P`, specify a measured observable or likelihood `p(y | gamma, eta, P)`, where `eta` includes calibration and nuisance parameters. For example a frequency encoding is `omega = kappa gamma`; `kappa` has units of inverse seconds and must be calibrated. Cosmic age, laboratory elapsed time, coherent interrogation time, dimensionless driving coordinates, and iteration count must remain distinct variables.

A useful pointwise criterion is

\[
\Pr_P(|\widehat\gamma_j-\gamma_j|\leq\delta_j)\geq1-q.
\]

Choose absolute tolerance, fractional tolerance, or tolerance relative to a measured local spacing **in advance**. The mean spacing `2 pi / log(T/(2 pi))` is not a lower bound on the smallest nearby gap. Demonstrating a high isolated zero also differs from resolving every zero in a contiguous prefix. If claiming all zeros through `T`, account for multiple comparisons or simultaneous coverage.

An accessible-height definition can be a supremum over a stated class of protocols meeting this criterion and a finite resource vector, with a separate definition for isolated targets versus a complete prefix. A nondecreasing resource-feasible set as cosmic age increases implies a nondecreasing *optimized* accessible height by set inclusion. That monotonicity is conditional and mathematical; whether usable physical resources are nested is additional physics. Accelerated expansion, loss of causal access, available free energy, decoherence, and finite calibration stability can invalidate simple extrapolation of local resources.

Treat model discrepancy and calibration bias explicitly. Adding squared biases and statistical variances is justified only with a specified risk definition or probabilistic model. It is not a substitute for identifying and controlling a systematic effect.

## Precision, resolution, and quantum bounds

The quantum Cramér–Rao framework bounds locally unbiased estimation variance by quantum Fisher information. For the illustrative pure-state Ramsey model `H_gamma = (hbar kappa gamma / 2) sigma_z`, with an initial equatorial qubit and interrogation duration `tau`, the information per independent trial is `F_Q = (kappa tau)^2`. Consequently

\[
\operatorname{Var}(\widehat\gamma)\geq \frac{1}{M\kappa^2\tau^2}
\]

for `M` independent trials under the regularity assumptions. This is a local lower bound on variance, not a guarantee that an estimator achieves it, not a resolution theorem for unknown numbers of lines, and not a global unwrapping protocol. Finite contrast, nuisance parameters, prior knowledge, decoherence, and calibration can change the inference. The example is our elementary application, not an asserted universal model of all zero experiments. [Braunstein and Caves, 1994](https://doi.org/10.1103/PhysRevLett.72.3439); [Giovannetti, Lloyd, and Maccone, 2006](https://arxiv.org/abs/quant-ph/0509179).

The familiar Fourier width `2 pi / tau` is a benchmark for specified windowing and reconstruction. Slepian and Pollak address time–band concentration, not a universal ban on estimating two nearby frequencies more precisely. [Slepian and Pollak, 1961](https://doi.org/10.1002/j.1538-7305.1961.tb03976.x).

Concrete quantum-sensing protocols can change the resolution behavior, including superresolution under an explicit noise model. This is a direct reason not to make the Rayleigh criterion a law of nature. [Gefen, Rotem, and Retzker, 2019](https://www.nature.com/articles/s41467-019-12817-y). External synchronization also changes the relevant available time resource; a sensor's coherence time need not alone set the spectral resolution. [Boss et al., 2017](https://arxiv.org/abs/1706.01754).

For an engineered scalar signal `g(E) = -zeta(1/2+iE)/(1/2+iE)`, a simple zero has derivative magnitude `|g'(gamma_j)| = |zeta'(rho_j)|/|rho_j|`. Linear error propagation would give `sigma_gamma ≈ sigma_g / |g'|` if this entire complex signal were measured with the assumed noise. A real-channel experiment needs its actual real-channel derivative. This illustrates why zero-dependent slope and signal normalization matter: a common fixed noise level does not imply one universal cutoff. It should not be reported as a refit of a particular experiment without its likelihood and calibration.

## What information and energy limits do—and do not—provide

Bekenstein's entropy–energy relation motivates finite information capacity for a specified bounded system under its assumptions. It does not directly identify an allowed largest prime or Riemann-zero ordinate. [Bekenstein, 1981](https://doi.org/10.1103/PhysRevD.23.287).

The Margolus–Levitin result constrains the rate of orthogonal state evolution using energy above the ground state. Turning this into a computation budget requires a defined operation model; turning it into a zero bound additionally requires an algorithmic cost lower bound. [Margolus and Levitin, 1998](https://arxiv.org/abs/quant-ph/9710043).

Lloyd's cosmological estimates quantify information-processing budgets under a cosmological model. They are useful motivation for resource accounting, but do not imply that the world has only a finite mathematical sequence of primes or only 4200 physical zeros. [Lloyd, 2002](https://arxiv.org/abs/quant-ph/0110141).

An explicit catalog of `J` independently tabulated values at `b` bits each needs roughly `J b` bits, whereas an algorithm plus an index can be much shorter. Therefore a storage bound is not generally an index bound: mathematical sequences are compressible. A hardware bandwidth constraint also depends on the encoding factor `kappa`; changing that factor trades range against sensitivity and does not reveal an invariant cosmic ordinate threshold.

## Prime reconstruction is a separate numerical task

The original formula for `pi(x)` omits essential prime-power/Möbius terms and offers an `O(x)` remainder too large to support the proposed precision. A cleaner exact starting point is the symmetrically summed explicit formula for the half-weighted Chebyshev function:

\[
\psi_0(x)=x-\sum_{\rho}\frac{x^\rho}{\rho}-\log(2\pi)-\frac12\log(1-x^{-2}),\qquad x>1.
\]

At prime powers the midpoint convention matters. Finite-height versions have truncation errors depending on `x`, distance to discontinuities, smoothing, and zero height. There is no general theorem that adding zeros gives monotonically, exponentially decreasing pointwise error or a unique universal sufficient cutoff. A reconstruction threshold at a chosen `x` and tolerance is not an experimental detectability bound. Avoid making a rigorous asymptotic remainder claim unless its hypotheses and primary mathematical reference are included in the actual manuscript.

The fine-structure variation amplitude is a different dimensionless observable from waveform noise, phase uncertainty, and relative prime-counting error. Setting all of them equal to `10^-4` requires a derived transfer relation. No such relation is supplied in the original argument. The arithmetic label “p-adic” likewise requires an actual p-adic construction; selecting primes in a multiplicative empirical scaling does not establish one.

## Update in the last six months

The targeted search did locate a directly relevant new experiment:

**Shijie Wei et al., “The Riemann Hypothesis manifested in dynamical quantum phase transitions,” Nature Communications 17, 8163, published 1 July 2026; DOI 10.1038/s41467-026-74935-8.** A five-qubit nuclear-spin proof of principle accompanies engineered quantum models. The large-index demonstrations through the trillionth zero are numerical simulations, not direct measurements of that many zeros. The paper explicitly rescales physical time by an independently specified energy scale. Its public source-data repository can support a concrete update while preserving this distinction. [Publisher article](https://www.nature.com/articles/s41467-026-74935-8); [author data and code](https://github.com/lqf2025/Riemann-data); [Zenodo record](https://doi.org/10.5281/zenodo.20590665).

Its preprint first appeared on 14 November 2025: [arXiv:2511.11199](https://arxiv.org/abs/2511.11199). Thus the journal publication is new in the search interval; the underlying idea is older. The original 2021 trapped-ion paper remains a separate experiment and should retain its actual date and correct authors: [He et al., npj Quantum Information 7, 109](https://doi.org/10.1038/s41534-021-00446-7).

For the revision, cite the 2026 paper's finite-system observations and resource-dependent simulations. Do not import its strongest interpretation as a proven physical origin of arithmetic or a certified computational speedup. In particular, normalized thermodynamic-limit signals require care: a normalization that diverges can make fixed-coordinate signals vanish generically, so finite-system minima and scaled observables should be distinguished from exact limiting zeros.

The search did not establish a new measured **cosmic-age dependence** or a universal experimental zero-height ceiling. Absence from this targeted search is not proof that no relevant work exists.

## Mathematical computation is a different evidence category

Platt and Trudgian rigorously verified RH for all positive ordinates through `3 × 10^12` using interval arithmetic. This is a mathematical computation on finite physical hardware, not a direct analog-spectroscopy experiment. It decisively precludes interpreting 4200 as a general bound on physically computable zero information, while leaving an explicitly restricted analog protocol open for study. No claim of a new record is necessary. [Platt and Trudgian, Bulletin of the London Mathematical Society 53, 792–797 (2021)](https://doi.org/10.1112/blms.12460); [author preprint](https://arxiv.org/abs/2004.09765).

## How the author's newly published dynamics paper fits

The non-autonomous-map paper can motivate a candidate aging law, but its iteration index is not established cosmic time. A response law from the map parameter to experimental noise, generator strength, or accessible resources is an additional assumption. The prior numerical transfer operator is not established as the self-adjoint Hilbert–Pólya Hamiltonian. Its result should be cited as numerical spectral correspondence, not strict isomorphism or proof of all local GUE statistics.

An honest revision can therefore discuss two nested levels: a reproducible protocol-limited observation model, and a separately labeled cosmic-aging extension. Demonstrating the former must not be presented as evidence for the latter. An inverse-log-squared aging ansatz may be retained if the observable it modifies, units, normalization, and free parameters are defined; its superiority or cosmic interpretation require independent data.

