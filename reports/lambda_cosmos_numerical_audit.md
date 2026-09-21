# Numerical and mathematical audit of the lambda-cosmos truncation update

Audit date: 2026-09-20. Source: `lambda-cosmos_v1_latex.pdf`, 22 pages,
SHA-256 `042eb4e99f790260af41f403c835e8a49977af61e69e9e1a119b05bb6f9661f1`.
The PDF preflight passed. Page numbers below refer to physical PDF pages.
This report concerns the zero-truncation material in Sections 3.3, 4.3 and S4;
it does not reanalyse the source's alpha or Hubble fits.

The newer manuscript improves some bookkeeping relative to the earlier
`riemann_zeros_Fixed.pdf`, but its universal physical-cutoff inference remains
unsupported. Its newer 200–400 million year scale deserves a more careful
treatment than the earlier, inconsistent 100 million year estimate: the new
scale can follow conditionally from specified logarithmic assumptions.

## 1. What the source actually adds

- Section 3.3 and Figure 4, pp. 10–11, explicitly distinguish the proposed
  laboratory height 73.3 from an interpolated zero index near 18.5. It still
  claims a physical collapse near index 80 without establishing one from the
  data. The current revision should not say this newer source simply calls
  height 73.3 “the 80th zero”; that would misdescribe the source.
- Section 3.3.2, p. 12, updates the old approximate cosmic height 4200 to 4397,
  gives a purported exact count of 3884 zeros, uses mean spacing 0.96, and
  estimates 200–400 million years to expose a new zero.
- Section 4.3, p. 14, proposes an observational look-back-time comparison using
  high-redshift absorption spectra. This is an additional research idea, not
  a measurement already supplied by the manuscript.
- Supplement S4, pp. 19–20, supplies a proposed phase-noise/Nyquist argument.
  The valid phase propagation and unsupported sampling/energy conversion
  steps need to be separated explicitly.

## 2. Counts and spacings: the new cosmic arithmetic is still not exact

The formula used by the source gives

\[
T_0=\frac{\ln(10^{60})}{\pi\sqrt{10^{-4}}}
   =4397.613593276566445\ldots.
\]

Calculations with `mpmath` 1.3.0 at 40 decimal digits give:

| Height | Leading smooth count, without 7/8 | Smooth count including 7/8 | Computed integer count N(T) |
|---|---:|---:|---:|
| 4397 exactly | 3884.474116 | 3885.349116 | 3885 |
| 4397.613593276566 | 3885.113850 | 3885.988850 | 3886 |
| 73.293559887943 | 16.991236 | 17.866236 | 18 |
| 76 exactly | 18.057255 | 18.932255 | 19 |

The source's “exactly 3884” comes from treating a smooth asymptotic count as
an exact integer, compounded by rounding the proposed height. Neither its
rounded nor unrounded height has that exact count. The integer results and
the bracketing zeros at the cosmic height were independently recomputed at
60 decimal digits. These are numerical checks, not formal interval proofs.

The neighboring ordinates are:

| Index j | Ordinate gamma_j |
|---|---:|
| 3883 | 4394.028511126826 |
| 3884 | 4395.060893844163 |
| 3885 | 4395.951097930616 |
| 3886 | 4397.223222687317 |
| 3887 | 4398.451719633070 |

Thus, for the unrounded proposed height, the next zero is index 3887 and
the remaining height is 0.838126356503, not a full mean gap and not a step
from index 3884 to 3885. At the rounded height 4397, the next zero is instead
index 3886; the rounding itself crosses one zero.

The quoted spacing is a reasonable **local mean-density approximation**:

\[
\overline{\Delta\gamma}(T_0)
 =\frac{2\pi}{\ln(T_0/2\pi)}=0.959127251332.
\]

It is not the exact gap between a selected neighboring pair. The smooth
count also has a constant term and an oscillatory remainder; rounding its
leading term cannot supply a universally correct integer count.
For the distinction between counting all strip zeros and computing their
critical-line ordinates, see [NIST DLMF, Section 25.10](https://dlmf.nist.gov/25.10).

## 3. The 200–400 million year scale is conditional, rather than an arithmetic impossibility

To differentiate the proposed height law, the source needs separate,
dimensionless definitions of its information scale and its noise variance.
One explicit realization, normalized to retain the source's present values,
is

\[
X(t)=10^{60}\frac{t}{t_0},\qquad
\epsilon(t)=10^{-4}
 \left[\frac{L_0}{\ln(t/t_*)}\right]^p,
\quad L_0=\ln(t_0/t_*).
\]

Here take the same illustrative values as the current revision:
`t0 = 13.8 Gyr`, `t* = 5.391247e-44 s`, and 31,557,600 seconds per year.
Then `ln(X0) = 138.155105580` and `L0 = 140.244226814`.
These two logarithms should not silently be identified: the actual age in
these Planck-time units is approximately `8.0778e60`, not `1e60`.

Writing `y = ln(t/t0)`, the normalized height is

\[
\frac{T(t)}{T_0}
 =\left(1+\frac{y}{\ln X_0}\right)
  \left(1+\frac{y}{L_0}\right)^{p/2},\qquad
\left.\frac{dT}{dt}\right|_0
 =\frac{T_0}{t_0}
  \left(\frac{1}{\ln X_0}+\frac{p}{2L_0}\right).
\]

If one instead sets `X=t/t*` exactly, the simpler law is
`T proportional to [ln(t/t*)]^(1+p/2)`, but then its present normalization
must also be adjusted consistently.

| Assumed noise evolution | Linearized wait for one mean-gap increase | Exact wait for that increase | Exact wait from T0 to gamma_3887 |
|---|---:|---:|---:|
| p=0: constant epsilon | 415.82 Myr | 422.15 Myr | 368.19 Myr |
| p=2: epsilon proportional to inverse log squared | 209.47 Myr | 211.06 Myr | 184.25 Myr |

Consequently, an order-of-magnitude 200–400 Myr result is compatible with
some versions of this phenomenological law. It is not a unique prediction
of it, and a mean-gap waiting scale is different from the wait to the next
specific ordinate. The source leaves the growth law for X, the definition
and normalization of epsilon, the physical meaning of T, and the
noise-to-resource mapping unmeasured. The apparent proximity to a
geological cycle is not independent validation of any of these assumptions.

There is also a convention issue when comparing this source with the new
paper: S4 treats epsilon as a **variance**, whereas the current revision's
delta_cos is an **ordinate-resolution scale**. A rule epsilon proportional
to inverse log squared cannot automatically be transferred to delta_cos
with the same exponent. Their relationship is itself a response-model
assumption. If epsilon were proportional to delta_cos squared, a
delta_cos inverse-log-square law would instead correspond to p=4 here.
That hypothetical identification gives a 140.00 Myr linearized mean-gap
scale (140.70 Myr using the exact height increment), rather than the p=2
value. It is only a convention comparison: equality between a variance in
logarithmic encoding coordinates and an ordinate-resolution variance has
not been established.

Recommended manuscript treatment: preserve this conditional derivative
calculation in an audit appendix or provenance discussion, but do not
restore an observational cosmic “next zero” forecast.

## 4. The look-back-time idea is not automatically a large effect

Under the explicit p=2 assumptions above, the proposed height changes only
logarithmically over astronomically large age ratios:

| Cosmic age | T(t)/T0 | Relative reduction from the present |
|---|---:|---:|
| 6 Gyr | 0.988068 | 1.19% |
| 3 Gyr | 0.978193 | 2.18% |
| 1 Gyr | 0.962643 | 3.74% |
| 0.5 Gyr | 0.952896 | 4.71% |
| 0.1 Gyr | 0.930455 | 6.95% |
| 0.001 Gyr | 0.867722 | 13.23% |

These are illustrative **age scenarios**, not new redshift observations or
a fitted cosmology. Converting them to redshift requires a specified
expansion history. In particular, the source's claim that an early-epoch
cutoff must be “drastically smaller” does not follow just from the specified
logarithmic schedule. A different choice of reference time or noise
evolution could give different effects, which makes these additional model
parameters rather than established predictions.

The proposed astronomical route also needs a forward model from atomic
transition spectra, environments and instrument response to an operational
zero-resolution observable. A qualitative similarity of spectral statistics
does not identify atomic transition energies with individual Riemann
ordinates. The current paper can retain the look-back-time comparison as a
conditional research direction, without claiming that quasar absorption
spectra already measure a zero truncation boundary.

## 5. S4: valid phase propagation, but no universal Nyquist derivation

Put `u=ln X`, and consider an explicitly chosen encoded component
`exp(i gamma u)`. If u has absolute standard deviation `sigma_u=sqrt(epsilon)`,
then at fixed gamma the phase standard deviation

\[
\sigma_\Phi=|\gamma|\,\sigma_u
\]

is valid. It does not, by itself, establish the following two subsequent
steps in S4.

**Sampling.** Uniform sampling of u at spacing Delta u has angular Nyquist
frequency `pi/Delta u`. The observation window length, approximately ln X
under a particular setup, controls spectral resolution; it does not specify
the sampling interval. S4 supplies no sampling grid or transfer function
that turns this theorem into `T sqrt(epsilon) <= ln X/pi`. That inequality
is an additional threshold choice. The bandlimit and sample-spacing
requirements are explicit in [Shannon's original paper](https://fab.cba.mit.edu/classes/S62.12/docs/Shannon_noise.pdf).

For example, a stated Gaussian encoding-jitter model instead gives

\[
\left|\mathbb E e^{i\gamma\delta u}\right|
 =\exp(-\gamma^2\sigma_u^2/2).
\]

A chosen visibility threshold `V_min` would yield the protocol-dependent
coherence scale `sqrt(-2 ln V_min)/sigma_u`. This derives neither the
specific factor ln X/pi nor a universal cliff. White noise, statistical
dephasing and aliasing are different mechanisms.

**Ordinate error.** S4.3 turns `sigma_Phi = 76*0.1 = 7.6` into an ordinate
uncertainty of the same size. These are different observables. For a simple
unwrapped phase estimator using the nominal u, the encoding-jitter
contribution would instead be

\[
\sigma_\gamma\simeq\frac{|\gamma|\,\sigma_u}{|u|}.
\]

With `u=ln(1e10)`, gamma=76 and sigma_u=0.1, this is 0.330064,
whereas the local mean ordinate gap is 2.520476. This illustration is not
a fitted model of the ion experiment and does not imply it achieves that
uncertainty; it exposes the missing conversion in S4. An assumption about
relative u error could give another formula, but would differ from S4.2's
stated absolute error. A quantum wavefunction overlap or energy-level
broadening conclusion additionally requires an actual Hamiltonian and
noise channel, which the explicit-formula phase alone does not supply.

Finally, the numerical substitution `X=1e60`, `epsilon=1e-4` is a chosen
normalization, not an independent measurement of a universal cosmic
precision. The relevant quantity cannot simultaneously be treated as an
alpha-drift amplitude, a variance in ln X, an absolute phase tolerance,
and an ordinate error without explicit physical conversion laws.

## Reproducibility and disposition

Machine-readable arithmetic is in
`results/lambda_cosmos_arithmetic.json`. It records numerical precision,
versions, constants, exact scenario definitions, neighboring ordinates and
the derived waiting and age scenarios. Counts used `mpmath.nzeros`; nearby
ordinates used `mpmath.zetazero`. The consistency checks support the
arithmetic above; they do not establish any physical cosmic bound.

Keep the newer source in the project's provenance and cite it as the
antecedent of the phase-noise and look-back-time proposals. Retain the
conditional ideas with the missing assumptions exposed. Correct the
integer count, qualify the mean spacing and waiting scale, and do not
reintroduce universal physical truncation, a confirmed noise wall, or a
connection to geological cycles as demonstrated results.
