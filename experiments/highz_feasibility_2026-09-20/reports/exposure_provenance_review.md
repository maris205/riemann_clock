# Individual ESPRESSO exposures: independent provenance and schema review

The 17 public FITS products are usable for conditional native-pixel consistency tests. They are **already extracted, wavelength-calibrated, sky-subtracted spectra**, rather than raw detector images or 17 independent gas-model determinations. The analysis must retain their native wavelength grids and transfer the published exposure-specific rejection masks before interpreting order/epoch differences.

## Sources and archive

The source data repository is [Murphy et al., ESPRESSO HE0515−4414](https://github.com/MTMurphy77/ESPRESSO_HE0515-4414), pinned by the download manifest to commit `09ee2cb32fb57ace9cda2943e997930bc7813366`. The associated [paper](https://arxiv.org/abs/2112.05819) describes DRS 2.2.3, 17 exposures, two traces per order, and additional clipping and atmospheric masks used during combination. Its final combined spectrum and the exposures contain the same photons.

The official [UVES_popler ESPRESSO reader](https://github.com/MTMurphy77/UVES_popler/blob/711b07badfc8209105d0de437f56f4cea1259563/UVES_r2Dspec_espresso.c), wavelength conversion, action replay, and atmospheric-mask routines were downloaded at commit `711b07badfc8209105d0de437f56f4cea1259563`. They are archived as `data/raw/exposures/provenance_UVES_*.c`, with SHA256 and URLs in `provenance_sources.json`. This is a current reader-source audit, not a claim that this exact later commit produced the historical combined spectrum.

The [ESO ESPRESSO DRS manual](https://ftp.eso.org/pub/dfs/pipelines/espresso/espdr-pipeline-manual-3.3.10.pdf) is also archived. It describes the same named S2D extensions but is version 3.3.10; the actual FITS headers establish that these observations used version 2.2.3. The [Dumusque et al. 2021 data-release appendix](https://arxiv.org/html/2009.01945v2) corroborates the extension conventions for the ESPRESSO 2.2.3 DRS adapted to HARPS-N. It is supporting documentation for the software format, not an additional ESPRESSO observation.

## Confirmed schema and operations

| Quantity | Finding | Required handling |
|---|---|---|
| Shape | Each extension is 170 × 9111; two adjacent rows represent the traces of each of 85 order entries. | Keep rows separate. Row `r // 2` is the order entry; `r % 2` is the trace. These are array indices, not physical echelle interference-order numbers. |
| Wavelength | `WAVEDATA_VAC_BARY` is vacuum wavelength in Å, already corrected to the barycentre. | No air conversion and no second BERV correction. |
| Flux | `SCIDATA` is extracted flux per pixel, not continuum-normalized. `S2D_SKYSUB_A` identifies the sky-subtracted science-fibre product. | Model the local continuum separately for each trace; do not subtract sky again. |
| Error | `ERRDATA` contains the flux standard uncertainty, not variance. The UVES reader uses it directly. | Divide residuals by this array, not its square root. A flux normalization must also divide the error by the same scale. |
| Pixel quality | `QUALDATA == 0` is the accepted convention. The UVES reader rejects every nonzero quality value and nonpositive errors. | Reject all nonzero flags plus nonfinite flux/error and nonpositive errors. Individual bit meanings were not decoded; this does not affect an all-nonzero exclusion rule. |
| DLL | `DLLDATA_VAC_BARY` gives wavelength width, not another flux correction. The official reader never reads this extension and does not divide the flux by DLL. | Midpoint boundaries from adjacent native wavelength values are valid here. A conversion to flux density would require changing flux, sigma and continuum consistently; it is unnecessary for the present local count-profile fit. |
| Blaze | Standard S2D products are deblazed; `S2D_BLAZE_*` denotes the separate non-deblazed category in the ESO DAS documentation. | Do not apply an additional blaze division. Remaining smooth throughput/continuum is fitted locally. |

The `helio = 0` and `vacwl = 0` settings in the UPL do **not** mean a second conversion should be applied to ESPRESSO. `UVES_wpol.c` explicitly returns the imported wavelengths for the ESPRESSO file type before the generic air/heliocentric conversions. Its generic source comment calls these spectra “heliocentric,” but the imported HDU is `WAVEDATA_VAC_BARY` and the code preserves it. All 17 UPL velocity-shift triples are zero, and `distort = 0`.

## Published masks: material correction found during review

The original UPL records 2,781 actions: 2,460 order clips, 304 order-continuum fits, 11 combined-continuum fits and six combined-spectrum clips. There are no unclip actions. Of these, **213 order clips overlap the three selected Fe II regions**: 30 at 2374, 118 at 2382 and 65 at 2600. The selected rectangles have flux limits far wider than any observed count, so they act as wavelength exclusions. The one-based UPL exposure and `ORDR` fields map to the exposure list and FITS row plus one. These relationships were independently verified against the original UPL and raw FITS rather than merely accepting the prepared mask JSON.

The original clips operate on redispersed trace pixels. Excluding native pixels whose wavelength bins overlap those published rectangles is a **conservative transfer**, not exact replay of the original redispersion and clipping. That transfer rejects 2,274 additional native samples over the selected exposures and line regions before atmospheric-mask exclusions. A nearest-pixel mask from the combined spectrum alone cannot reproduce exposure-specific clipping.

The supplied atmospheric mask also needs exposure-specific treatment. Shift its observer-frame wavelength intervals to the frame of the exposure using the stored BERV before checking native-bin overlap. The final analysis expands the transferred boundaries by ±0.25 km/s and excludes 433 additional samples after the manual-mask transfer. In the present three regions, the atmospheric intervals overlap the blue edge of Fe II 2600 for the nine 2018 exposures and the 2019 exposure; there are no intersections at 2374 or 2382. For example, the first exposure has an interval ending near 5590.9168 Å, inside the fitted range that starts at 5590.810 Å. A combined-spectrum pixel may remain valid when other exposures provide unmasked flux, so the combined mask alone is insufficient. The original UVES atmospheric-mask routine uses its computed heliocentric velocity, whereas using the exposure BERV is the direct native-barycentric transfer; this is another reason not to call the adaptation exact UPL reproduction.

## Interpretation limits

The shared gas template was learned from the same combined observations. Therefore exposure, trace and order comparisons are conditional consistency checks. Their absolute means are neither independent confirmations of the combined-spectrum fit nor measurements of temporal evolution of physical constants. Changing from a full simultaneous gas re-fit to a fixed-template estimator changes the estimator; compare the native exposures with a combined-spectrum result obtained by the **same fixed-template method**.

Different exposures have disjoint photon samples, but share calibration, extraction, atomic inputs and the data-derived template. The working independent-error likelihood does not establish that extraction/sky-subtraction noise and calibration errors have zero covariance. Native sampling avoids adding the final combination's interpolation correlation; it does not prove that all original pixel errors are independent.
