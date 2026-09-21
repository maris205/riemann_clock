# Independent review of the UVES SQUAD archive expansion

Reviewer scope: source-based quality policy, catalogue bookkeeping, and independent checks of the forthcoming extraction products. This review does not certify a cosmological measurement or a new absorption-model fit. Initial review: 2026-09-21.

## Source and bookkeeping checks

The local SQUAD DR1 inputs are pinned to commit `a0cdc8e7b99f2b01a45d919988af9d60e6447d19` of the [authors' data repository](https://github.com/MTMurphy77/UVES_SQUAD_DR1/tree/a0cdc8e7b99f2b01a45d919988af9d60e6447d19). I read the original `Paper/src/paper_accepted_2018-10-16.tex`, `Notes_FITS_Files.txt`, master catalogue, DLA catalogue, and the local NIST precision-wavelength extract. Independent CSV parsing gives **475 master-catalogue quasar rows**, **467 quasars with final spectra**, **132 DLA sightlines**, and **155 individual DLA absorbers**. Multiple absorbers on the same sightline share instrumental calibration and must not be counted as independent calibration experiments.

The 155 absorbers are an identified, targeted DLA sample. They are not a complete blind sample of all Fe II absorbers in 467 spectra. A screen of these 155 systems must not be described as having exhausted the archive. The earlier HE0515 and Q0347 examples should be explicitly excluded from the new-target count.

## Recommended frozen screening policy

1. Expand comma-separated absorber redshifts to one row per absorber, preserving its sightline identifier and catalogue source. Use a declared, consistent precision vacuum wavelength basis; integer line labels are not wavelengths. The implemented nine-line screen uses the archived Murphy–Berengut composite entries in `MM_VPFIT_2013-11-10_noiso.dat`: Fe II 1608, 1611, 2249, 2260, 2344, 2374, 2382, 2586, and 2600. I independently checked all nine wavelengths and oscillator strengths against that source. The weaker 1611/2249/2260 lines can be useful diagnostics but do not count toward the frozen requirement of three lines with f ≥ 0.01. They share the ground level, but neither this fact nor a profile-centroid comparison supplies the unknown clock response coefficient.
2. Require the entire fixed window, initially ±250 km/s in logarithmic velocity, to lie within one catalogue coverage interval. Do not use only WavStart/WavEnd: UVES spectra have internal gaps. Catalogue coverage itself smooths over some short invalid runs, so this is only a pre-screen.
3. Require the full window to be redward of the quasar Lyα forest. A fixed 3000 km/s redward margin follows the source paper's conservative survey convention. Separately flag systems within 3000 km/s of the quasar emission redshift. Quasar emission redshifts are approximate and this is a screening convention, not an exact physical boundary.
4. Record explicit telluric masks and apply them to the entire window. The source paper identifies O2 near 6300, 6880, and 7620 Å, H2O bands around 6470–6600, 6830–7450, 7820–8620, and 8780–10000 Å, and ozone below roughly 3300 Å. Broad exclusion intervals are conservative policy choices and must not be represented as precise transmission measurements. Account for the fact that tellurics move with heliocentric correction, e.g. by a declared fixed velocity padding. A line outside these masks is not guaranteed telluric-free.
5. Flag sky residuals near 5578 and 5891/5897 Å and recurrent DR1 artefacts near 4716, 4744, 5240, and 5580–5800 Å. The latter are not present in every target, so they may be retained for inspection rather than automatically rejected. Report their flags.
6. Use the nearest catalogue CNR entry only as a proxy. The five values refer to 3500, 4500, 5500, 6500, and 7500 Å and are normalized to **2.5 km/s per pixel**, while spectra may have other dispersions. A proxy threshold of 20 cannot certify the noise at a particular line, especially beyond the tabulated wavelength range. Resolve this using actual spectrum pixels or the fine CNR map.
7. Rank and select candidates using declared quality quantities, before inspecting any apparent shift or model-preference result. A download budget such as the top five is an implementation limit, not evidence that the sixth candidate is unusable. Report all pre-screened targets and reasons for rejection, including ties.

## Required checks after download

- Verify source URLs, byte counts, SHA-256 hashes, FITS identity, and the expected nine-array semantics. Reconstruct the vacuum–heliocentric log wavelength from `CRVAL1`, `CRPIX1`, and `CD1_1` with the FITS one-based pixel convention. Keep error and expected-fluctuation arrays distinct; they are not interchangeable definitions of covariance.
- Validate every window using finite flux, positive errors, positive expected fluctuation where used, and pixel status 1. Record expected/available pixel counts, valid fraction, longest invalid run, and contributing-exposure counts. A catalogue interval can conceal up to 10 km/s gaps under the archive's own map definition.
- Estimate local continuum CNR in predeclared sidebands, reporting missing or absorbed sidebands. Check automatic-continuum bias and retain continuum freedom in any later physical fit. The source paper explicitly reports underestimated low-flux uncertainty and remaining narrow cosmic-ray artefacts; saturated cores need separate attention.
- Record a weak-line detection diagnostic and saturation fraction. Three covered windows need not mean three detected Fe II lines. Avoid filtering individual pixels by their observed residual sign, and do not describe automatic screening as definitive blend identification.
- Confirm the same velocity structure is plausibly present in at least three usable transitions, including a weaker line if available. Strong and weak lines can have different apparent flux centroids under perfectly constant atomic frequencies because of saturation and unresolved component weighting. A common physical component model is necessary before interpreting relative line offsets.
- Treat the nominal resolving power as a slit-illumination estimate. The source paper warns that the actual quasar resolving power can be larger; do not fix it to the catalogue value and interpret the resulting profile mismatch as new physics.

## Precision and inference limits

DLA redshifts are often quoted to three decimal places. At redshift z, rounding by 0.0005 corresponds to approximately `c*0.0005/(1+z)` in velocity, tens of km/s. These redshifts locate search windows, but cannot define a precision velocity zero. Fit a common absorber redshift or velocity nuisance; catalogue-centre offsets are not drift measurements.

The NIST extract provides identification and laboratory references. Its Fe II wavelength uncertainties correspond to roughly 13–15 m/s for the selected lines, and may share calibration or level-energy covariance. They do not encode the same isotope-composite convention as every published many-multiplet atomic list. Mixing conventions can create artificial offsets. Any sub-km/s analysis must document the chosen consistent isotope and wavelength basis.

The present archive expansion can establish that usable observations already exist, identify targets for component modelling, and quantify data limitations. It cannot itself establish alpha variation, a cosmic frequency drift, or a preference for `1/ln²(t/t*)`. Such conclusions additionally require a specified physical response, treatment of instrumental distortions and atomic-data covariance, and comparisons of time laws under the same sample-selection and nuisance model.

## Output validation status

The independent script `code/archive_quality_validate.py` reconstructs the full screen from the two original CSV files and the independent no-isotope atomic file, without importing the implementation under review. It verifies all 1395 line windows (155 absorbers × nine transitions), their laboratory values, observed wavelengths, full-window coverage and gates, catalogue CNR proxies, per-absorber gates, ranks, and selected sightlines. At the first frozen-screen snapshot, all **5899 mechanical checks passed**. The machine-readable result is `results/archive_expansion/independent_catalogue_validation.json`.

The screen yields **36 metadata-eligible absorbers**, with the following first five distinct sightlines, in the reproduced ranking:

| Sightline | Catalogue absorber redshift | Third-highest strong-line CNR proxy |
|---|---:|---:|
| J133335+164903 | 1.776 | 135 |
| J220852-194359 | 2.076 | 72 |
| J235129-142756 | 2.279 | 64 |
| J040718-441013 | 2.595 | 63 |
| J004131-493611 | 2.248 | 62 |

The source uses broader H2O ranges than the major-band exclusion intervals frozen in the implementation. I recommended retaining those ranges as explicit warnings, rather than relabelling every eligible line as atmosphere-free. The implementation has agreed to add them without making a shift-dependent selection change. This review will rerun after the final metadata edits.

The expanded scope downloaded **all 32 distinct sightlines containing the 36 metadata-eligible absorbers**. I independently opened the original archived FITS members, checked their SHA-256 hashes and byte counts, reconstructed the log wavelength on an algebraically independent basis, and compared all nine arrays, validity masks, exposure metadata, dates, and endpoint wavelengths with the processed products. I also independently reconstructed the secondary broad-H2O screen over all 155 absorbers. It leaves **seven metadata-eligible absorbers**, and its three additional download choices are correctly reproduced. These checks and the complete scope/source manifests pass **939 mechanical checks** in `results/archive_expansion/independent_data_validation.json` via `code/archive_data_quality_validate.py`.

The secondary screen adds a 30 km/s padding to the broader source-paper water bands. It preserves the original five targets and records the amendment before the line-profile diagnostic calculations. This is a transparent screen amendment, not a preregistered blind cosmological experiment. The relevant source subsection is **5.2.2, Telluric features**.

I also inspected the explicit wavelength-shift/slope settings in all 32 downloaded UPL files with `code/archive_upl_audit.py`. **None has a nonzero recorded VSHT shift or slope.** This is a useful archive-processing flag, not a complete independent calibration audit: a zero record does not prove every calibration operation absent. The source paper's subsection 5.2.9 documents instrumental distortions at scales relevant to a 100 m/s test, and states that only a minority of spectra receive additional corrections. These candidate coadds therefore require dedicated calibration assessment before supporting such a precision inference.

The completed pixel diagnostics were independently reconstructed directly from the 32 original FITS files, without importing `archive_quality.py`. The review covers every scalar equivalent width and diagonal error, actual pixel-edge integration, coverage and invalid runs, local inverse-error proxy, sideband flux, low-flux and moderately absorbed pixels, AOD normalization/floor masks and integrals, saved pixel indices and window arrays, atmospheric/artifact flags, table resolving powers, and the full readiness gates. **All 11,529 mechanical checks pass** across **324 windows**. Reproducer: `code/archive_diagnostics_validate.py`; detailed output: `results/archive_expansion/independent_diagnostics_validation.json`.

Independent reconstruction found a small integration normalization error during review: the nominal `UP_DISP` is 2.5 or 1.3 km/s, whereas the actual logarithmic grid spacing is `c ln(10) CD1_1`, lower by approximately 4.17 or 2.17 parts per million. The implementation now uses the latter for exact pixel widths and AOD integration, and records it separately. This correction changes no detection or readiness decisions.

The shared-absorption test was also clarified: it evaluates interpolated per-pixel absorption/error values at 5 km/s grid locations, not independent bin-averaged significances. The final implementation requires both native bracketing pixels to be valid and documents the requirement of five moderately absorbed pixels in a relatively weak transition. I independently confirmed that excluding interpolation across invalid runs changes no shared-grid counts or final decisions in this sample. The thresholds remain **detection proxies conditional on diagonal errors and the supplied continuum**, not calibrated probability statements.

Six absorbers pass the computational readiness gate:

| Sightline | z_abs | Detected strong transitions passing the strict screen | Shared grid locations | Edge-absorption warnings |
|---|---:|---|---:|---|
| J004131-493611 | 2.248 | 1608, 2374, 2382 | 8 | none |
| J053007-250329 | 2.141 | 1608, 2374, 2382 | 29 | none |
| J064326-504112 | 2.659 | 1608, 2374, 2382 | 21 | 1608 |
| J225719-100104 | 1.836 | 2344, 2374, 2382 | 25 | all three |
| J232128-105122 | 1.629 | 2344, 2374, 2382, 2586 | 7 | none |
| J233156-090802 | 2.143 | 1608, 2374, 2382 | 29 | all three |

These six are **candidates for developing a physical component model**. They are not six established constant-variation measurements. In particular, edge-absorption warnings indicate a need to revisit the rounded catalogue centre and finite fitting interval; broad atmospheric avoidance does not certify freedom from narrow tellurics or blends; and all require instrument calibration assessment. The strict-screen seventh candidate fails actual quality requirements, which demonstrates why catalogue selection alone was insufficient.

I visually inspected the J232128−105122 profile grid: the four strong selected transitions have a plausible shared dominant absorption complex and the weaker 2260 transition is visible. The stronger lines have nearly black cores, while 2374 retains flux. Additional weaker structure, especially in 2344, must be modelled before precision offsets can be compared. This supports its status as a modelling candidate without making any frequency-variation claim.

**Final review disposition:** implementation and reproducibility checks pass for the stated archive-readiness task. The numerical checks do not validate a physical drift, atomic response, fine-structure-constant change, or logarithmic cosmic-time law. The report's boundary between metadata screening, actual absorption detection, and future physical modelling must remain explicit.


## Final narrative and figure review

I reviewed `reports/archive_expansion_results_cn.md`, its generator `code/archive_report.py`, the rendered `archive_line_readiness.png`, and the J004131−493611 profile grid after the completed numerical validation. No blocking scientific wording issue was found. The reported 155 catalogue absorbers, 36 preliminary passers, 32 downloaded sightlines, 324 actual windows, seven strict catalogue candidates, and six computational candidates agree with the verified products. The six candidate rows, relatively weak-line requirements, and edge-absorption warnings match the numerical results.

The report explicitly limits its scope to the known DLA catalogue and does not imply a complete search of all metal absorbers in 467 spectra. It distinguishes historical photon observations from new analysis, treats VSHT zero records cautiously, and does not claim calibrated frequency changes or a cosmic-time fit. Its reference to this validation report does not turn the 18,367 mechanical checks into independent scientific evidence.

The rendered heatmap has 36 rows and nine transition columns; its six starred systems match the candidate list. The numeric cell labels correctly refer to median inverse normalized error per native pixel. Its four colours separate metadata rejection, failed native quality, broad-water warning, and detected lines outside broad water bands. The footnote explicitly preserves calibration, atmosphere, blending, saturation, and LSF limitations. Labels and legend are legible.

Two nonblocking editorial suggestions were sent to the report author: prefer “32 条不同视线” over “32 个独立视线” to avoid suggesting statistical calibration independence, and replace the future wildcard profile-report reference with the concrete pilot report when delivered. Physical-pilot fit validation belongs to that task's separate reviewer.
