# High-redshift spectroscopy: public-data pilot and conditional experiment design

**Actual public spectra are available and analyzed.** Beyond the initial two-sightline pilot, this package includes all 17 public ESPRESSO S2D exposures and 32 additional UVES SQUAD spectra / 36 screened absorbers. Successive physical-model and instrument controls are recorded by analysis stage. Clearly labelled synthetic injections and sample-size calculations remain separate. It does not measure cosmic evolution or a Riemann-zero cutoff.

## Focused spectroscopy manuscript and current follow-up

The [new single-column spectroscopy companion](../../paper_spectroscopy/feii_relative_frequencies.pdf) presents the astronomy methods and results independently of the [preserved 35-page Riemann-zero and cosmic-time manuscript](../../paper/riemann_clock.pdf). Its [source](../../paper_spectroscopy/main.tex) combines model sensitivity, native-order response controls, a common 2382/2374 descriptor and cross-age identifiability. Build it from this experiment directory with:

```bash
python ../../code/build_spectroscopy_paper.py
```

The command uses saved outputs and leaves the original manuscript unchanged. Full nonlinear analyses can be substantially slower than rebuilding the paper; individual reports document their inputs, frozen choices and validation scope.

- [Current follow-up summary](reports/publication_followup_results_cn.md)
- [Centred asymmetric order-response models](reports/order_response_results_cn.md) and [residual/trace diagnostics](reports/order_response_residuals_cn.md)
- [Common-pair conditional profiles and newly located local solutions](reports/common_pair_results_cn.md)
- [Cross-age design and identifiability](reports/cross_age_identifiability_cn.md)
- [Historical calibration products and remaining gaps](reports/calibration_readiness_cn.md)

The specified centred response family retains a 2600 order contrast near −59 m/s; this does not identify its unique cause. Profile searches expose local-solution dependence, so the earlier campaign values below remain historical comparisons. Cross-age calculations with assumed amplitudes are design tests, not measured cosmic evolution. The original hypothesis manuscript and synthetic forecasts retain their separate scope.

## Earlier completed campaign

See the [complete Chinese results](reports/available_data_campaign_results_cn.md) for the combined covariance/core re-fits, fixed gas-architecture sensitivity, individual exposures and orders, archive screen, two new physical fits and UVES continuation audit. Five-shift improvements under the joint masks are 23.69 and 6.70; the wider mask retains only 16% of the local information for the specified coherent strong-line direction. The same Fe II 2600 transition differs by about 60 m/s across adjacent orders under the conditional exposure estimator. Both new absorber pilots admit adequate conventional fits; the recommended-error controls give improvements of 6.16 and 5.50 for four and two shifts. These results do not establish physical time variation.

- [Combined-control validation](reports/joint_controls_independent_review.md) and [mask information](reports/mask_information_cn.md)
- [Gas architecture](reports/gas_structure_results_cn.md)
- [Exposure results](reports/exposure_results_cn.md) and [figure](results/exposures/exposure_consistency.png)
- [Archive results](reports/archive_expansion_results_cn.md) and [quality heatmap](results/archive_expansion/archive_line_readiness.png)
- [New physical fits](reports/archive_profile_results_cn.md) and [recommended-error control](reports/archive_profile_expected_noise_cn.md)
- [UVES strict continuation](reports/uves_completion_results_cn.md)

Reproduction commands for the new branches are provided at the end of this file. Earlier sections describe preserved stages of the investigation.

- [Full experiment protocol in Chinese](experiment_protocol_cn.md): physical assumptions, samples, controls, blind analysis, power, stopping rules and a current-to-2036 roadmap.
- [Data-readiness audit](reports/public_spectra_audit.md) and [actual observed spectral windows](figures/public_feii_windows.png).
- [Atomic response audit](reports/atomic_response_audit.md) and [instrument roadmap](reports/instrument_roadmap.md).
- [Conditional sensitivity figure](figures/conditional_sensitivity.png), [synthetic pixel-level pilot](figures/pixel_noise_pilot.png), and [machine-readable forecast](results/forecast_summary.json).
- [Independent forecast validation](reports/forecast_validation.md) and [pixel-pilot validation](reports/pixel_pilot_validation.md).
- [Independent protocol review](reports/protocol_review.md).
- [New actual-data results and outlook, Chinese](reports/existing_data_results_and_outlook_cn.md).
- [Observed measurements](results/observed_profile_metrics.csv), [registration configurations](results/observed_registration.csv), and [data-informed design](results/data_informed_design.json).
- [Observed-metrics independent review](reports/observed_metrics_independent_review.md), [registration review](reports/registration_validation.md), and [manuscript review](reports/highz_manuscript_review.md).

The initial two downloaded sightlines were HE 0515−4414 (reference absorber z=1.1508) and Q0347−383 (z=3.025). The historical reduced, coadded observations are not raw photon records or new 2026 exposures. The high-redshift spectrum has poor signal-to-noise in its two longest-wavelength Fe II windows, so it does not supply the six equally precise lines used in the separate planning forecast.

The central unresolved requirement is an **independently specified physical response** connecting the cosmic-time hypothesis to atomic levels or line shapes. The forecast's arbitrary response coefficients are only method-development inputs. A deterministic centroid-shift response is an extra assumption, not a derivation from the manuscript's uncertainty component. Alpha measurements, generic broadening, and astronomical line-spacing histograms cannot be renamed Riemann-cutoff measurements.

## Reproduce

The [public Hugging Face mirror](https://huggingface.co/datasets/dnagpt/riemann-clock-spectra)
provides the 52 downloaded spectrum products and 424 processed files. From this
experiment directory, `python ../../code/fetch_huggingface_spectra.py` restores the
large spectrum files at a pinned revision, verifies their checksums and preserves
the supplied manifests. It requires no authentication or additional Python packages.

Use Python 3.12 and the packages in [requirements.txt](requirements.txt), preferably in a virtual environment. Processed spectra, saved fits, small source snapshots and source manifests are included. Large third-party FITS files and spectrum archives are omitted from Git; restore them using the [repository contents and recovery guide](../../reports/repository_contents.md) before raw-data checks or new fits. Fetch scripts retrieve missing files and do not execute upstream code. Some retrieval scripts also regenerate manifests or processed arrays, so preserve the supplied snapshots before rerunning them.

```bash
python code/fetch_public_spectra.py
python code/fetch_atomic_data.py
python code/plot_public_spectra.py
python code/pixel_noise_pilot.py
python code/forecast_sensitivity.py
python code/validate_forecast.py
python code/measure_observed_profiles.py
python code/register_observed_profiles.py
python code/design_from_observed_data.py
python code/prepare_highz_manuscript.py
python ../../code/build_paper.py
```

The independently written pixel-pilot review is a saved audit, not an additional standalone test script. Source hashes are in [the spectral manifest](data/source_manifest.json) and [the atomic manifest](data/atomic/source_manifest.json). Numerical output retains seeds, assumptions and units. Scientific data and figure licenses retain their original notices; see the audits before redistribution.

The actual-profile analysis, two additional figures, a measurement table and an observing outlook are now included in the manuscript under `../../paper`. The complete experimental protocol remains here as supporting material. All new files remain within `riemann_clock`.

## ESPRESSO conventional null test

The added analysis now performs an actual physical baseline comparison on the public ESPRESSO HE0515−4414 spectrum: six Fe II transitions, 2,931 retained pixels, 45 shared gas components, isotope-resolved Voigt absorption, instrumental convolution, and jointly re-fitted gas and continuum parameters. This is a historical 2018–2020 coadd, not a new 2026 observation or an independent absorber.

The restricted conventional null is improved by five relative line shifts: **Δχ² = 51.71** in the selected primary comparison. Strong lines 2382 and 2600 favour offsets around −0.12 and −0.11 km/s relative to 2374. Several controls retain this conditional preference; removing those two strong lines reduces the improvement to 7.09 for three added shifts. The result is a residual consistency issue under specified assumptions, **not an established variation beyond conventional explanations, an alpha measurement, or an inverse-log-squared time-law detection**. The additional opacity/zero, strong-line removal, wider-bound and velocity-split checks are exploratory follow-ups.

- [Full Chinese results and interpretation](reports/conventional_null_results_cn.md)
- [Physical-model figure](figures/espresso_conventional_null.png) and [robustness comparisons](figures/espresso_null_robustness.png)
- [Selected comparison table](results/espresso_null/comparison.csv) and [initialization provenance](reports/espresso_null_run_provenance.json)
- [Source/mask/convention audit](reports/espresso_null_method_audit.md)
- [Independent numerical validation](reports/espresso_null_numerical_review.md), [saved-fit replay](results/espresso_fit_results_validation.json), [local injection controls](reports/espresso_null_local_controls_cn.md), and [velocity-split diagnostics](reports/espresso_velocity_split_local_controls_cn.md)

All third-party sources are pinned and hashed. Original budget-limited attempts are retained; the comparison selects converged continued fits. A local convergence flag does not certify the global optimum. Adjacent-pixel covariance and full model-selection uncertainty remain unmeasured. Local tangent Monte Carlo and injection recovery are explicitly separate from the constrained nonlinear comparison.

From this directory, a full replay uses the existing dependencies and cached inputs. The nonlinear fits take substantial CPU time; the original results are already included. Avoid running multiple copies with the same output names.

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python code/fetch_espresso_null_sources.py
python code/espresso_conventional_null.py --name null --max-nfev 350
python code/espresso_conventional_null.py --name alternative_from_published --free --max-nfev 350
python code/espresso_conventional_null.py --name alternative --free --start results/espresso_null/null.json --max-nfev 350
python code/espresso_conventional_null.py --name null_cross --start results/espresso_null/alternative_from_published.json --max-nfev 350
python code/espresso_null_robustness.py
python code/espresso_null_followup_controls.py
python code/polish_espresso_null_fits.py
python code/validate_espresso_null.py
python code/validate_espresso_fit_results.py
python code/espresso_null_local_controls.py --baseline null_cross
python code/espresso_velocity_split_controls.py --baseline null_cross
python code/prepare_espresso_manuscript.py
python code/write_espresso_null_report.py
python ../../code/build_paper.py
```

The export script stops if a selected comparison still lacks optimizer convergence. The raw and continued fit files, objective values, bounds and local-degeneracy diagnostics allow independent scrutiny rather than hiding failed attempts.

## Continuum covariance and frozen-core follow-up

The completed follow-up is documented in [the Chinese report](reports/conventional_followup_results_cn.md) and [the UVES report](reports/cross_instrument_results_cn.md). New continuum controls use 25 disjoint intervals / 13,750 pixels. The actual nonlinear GLS improvement is 29.24 for five additional shifts; separate frozen core masks give 41.66 and 13.26. These do not establish physical frequency changes or a cosmic-time law.

From this experiment directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/empirical_noise_controls.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/espresso_saturation_controls.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/polish_espresso_followup.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/espresso_empirical_gls.py
python code/prepare_followup_manuscript.py
python code/write_followup_report.py
```

The last two commands rebuild this earlier stage's manuscript inputs, figure and report from saved results. The polishing script retains failed attempts and reuses recorded continuations if present. The next stage now includes combined core/covariance controls, one alternative gas architecture and individual-exposure checks. Full extraction covariance and independently calibrated instrumental response remain unresolved.

## Reproduce the available-data extension

The saved fit arrays, processed inputs and source manifests are included; large original FITS files and spectrum archives must first be restored using the [recovery guide](../../reports/repository_contents.md) for checks that read or hash them. Source manifests pin the ESPRESSO and SQUAD repositories, preserve download URLs and hash each retrieved file. New numerical fits can take substantial CPU time. Run only one copy per output directory and preserve the supplied snapshots before re-running fit-producing scripts. No upstream downloaded analysis code is executed by these commands.

From this experiment directory, rebuild reports and manuscript inputs from saved results:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python code/prepare_available_data_manuscript.py
python code/exposure_report.py
python code/archive_report.py
python code/archive_profile_report.py
python code/uves_completion_report.py
python ../../code/build_paper.py
```

Independent saved-output checks, each preserving optimizer-status limitations:

```bash
python code/joint_controls_validate.py
python code/gas_structure_validation.py --saved-only
python code/exposure_validate.py
python code/archive_quality_validate.py
python code/archive_data_quality_validate.py
python code/archive_diagnostics_validate.py
python code/archive_profile_validation.py
python code/archive_profile_expected_validation.py
python code/uves_completion_validate.py
python ../../code/verify_available_data_delivery.py
```

The fitting entry points are `code/joint_controls_fit.py`, `code/gas_structure_refit.py`, `code/exposure_analysis.py`, `code/archive_profile_pilot.py`, `code/archive_profile_refine.py`, and `code/archive_profile_expected_noise.py`. Catalogue retrieval and quality processing use `code/archive_screen.py`, `code/archive_expand_all.py`, `code/archive_strict_screen.py`, and `code/archive_quality.py`; the order and fixed rules are documented in the archive report. Exact UVES continuation arguments and preserved starts are in its scope/results reports. The 40-component design and source hashes are saved before fitting. `code/mask_information.py` computes the separate fixed-parameter local information diagnostic, not another nonlinear fit.
