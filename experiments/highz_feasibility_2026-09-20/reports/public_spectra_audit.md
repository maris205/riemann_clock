# Public high-redshift spectra: actual data-readiness audit

Audit date: 2026-09-20. This is a reproducible archive pilot, not a measurement of cosmic evolution or Riemann-zero truncation.

## What is already available

**Yes: useful public, calibrated spectra exist now.** This pilot actually downloaded and decoded two UVES SQUAD DR1 spectra, plus the source catalogue and format documentation, rather than stopping at archive links. Downloaded scientific/reference sources total approximately 24.9 MB. The 475-row master catalogue contains 467 entries with final spectra (positive `NumExp`), consistent with the release. The catalogue counts quasars, not 467 independent absorbers or 467 suitable precision tests. [ESO release](https://www.eso.org/sci/publications/announcements/sciann17192.html), [author data repository](https://github.com/MTMurphy77/UVES_SQUAD_DR1), [Murphy et al. 2019](https://doi.org/10.1093/mnras/sty2834).

| Downloaded target | Quasar emission redshift | Reference absorber redshift | Observation dates in FITS exposure metadata | Exposures / exposure hours | Valid decoded pixels |
|---|---:|---:|---|---:|---:|
| HE 0515−4414 / J051707−441055 | 1.710 | 1.1508 | 1999-12-14 through 2009-02-11 | 51 / 50.46 | 279787 of 283993 |
| Q0347−383 / J034943−381030 | 3.205 | 3.025 | 1999-12-16 through 2009-09-23 | 29 / 39.28 | 241828 of 268238 |

The absorber redshifts come from published studies / the DLA catalogue; they are **not** inferred by equating the quasar redshift with the absorber redshift. HE 0515−4414 is an established metal-absorption precision target; Q0347−3819 has independently documented molecular and metal absorption at about 3.025. [Kotuš et al.](https://doi.org/10.1093/mnras/stw2543), [Levshakov et al.](https://doi.org/10.1086/324722).

These are historical coadded spectra. File retrieval in 2026 is not acquisition in 2026. The exposure table lists detector arms separately, so its 138 and 87 rows must not be mistaken for 138 and 87 independent observing exposures. Counts above use the FITS `UP_NEXP` values, consistent with the master catalogue.

## Actual arrays and processing conventions

Files are under `data/raw/` (unaltered author tar archives) and `data/processed/` (decoded NPZ arrays and metadata). The final spectra are coadded reduced data, **not raw photon records or unreduced detector images**.

Both spectra use vacuum, heliocentric wavelengths in angstroms and a logarithmic grid with native dispersion 1.3 km/s per pixel. Wavelengths were reconstructed from `CRVAL1`, `CRPIX1`, and `CD1_1`; both endpoints were independently checked against `UP_WLSRT` and `UP_WLEND`. The arrays have nine source channels: continuum-normalized flux, normalized statistical error, normalized expected fluctuation, continuum in arbitrary units, source status, contributing pixels before and after clipping, and chi-squared diagnostics before and after clipping. No flux fit, continuum refit, wavelength refit, or resampling was performed here.

`valid` means source `status==1`, finite flux/error, and positive error. It does **not** mean free of telluric absorption, sky residuals, blends, saturation, continuum error or wavelength-calibration distortion. Original redispersion creates correlations between nearby pixels; a full pixel covariance matrix is not supplied. The authors distinguish their statistical-error array from their expected-fluctuation array, recommending attention to the latter when fitting profiles. Precision work needs exposure-level checks, realistic covariance and calibration controls. [Pinned FITS documentation](https://github.com/MTMurphy77/UVES_SQUAD_DR1/blob/a0cdc8e7b99f2b01a45d919988af9d60e6447d19/Notes_FITS_Files.txt).

Decoded full-spectrum files:

- `data/processed/J051707-441055_squad_dr1.npz`
- `data/processed/J034943-381030_squad_dr1.npz`

NPZ keys: `wavelength_vacuum_heliocentric_AA`, `flux_normalized`, `error_normalized`, `expected_fluctuation_normalized`, `continuum_arbitrary_units`, `status`, `valid`, `contributing_pixels_before_clip`, `contributing_pixels_after_clip`, `chi2_before_clip`, `chi2_after_clip`. The corresponding metadata JSONs preserve all FITS header and table information. Individual Fe II-window NPZs also include `velocity_km_s`.

## Six Fe II windows: a concrete feasibility check

The twelve **real observed** windows are displayed in [the PNG figure](../figures/public_feii_windows.png) and [vector PDF](../figures/public_feii_windows.pdf), produced by `code/plot_public_spectra.py`. These panels contain no model fit or injected signal; the two poorest windows retain all observed values and error bands on explicitly expanded y-axes.

We checked approximate Fe II rest wavelengths 1608.4509, 2344.2128, 2374.4601, 2382.7642, 2586.6493 and 2600.1725 Å, mapped by the published absorber redshift. Windows are ±200 km/s for HE 0515−4414 and ±150 km/s for Q0347−3819. These windows do not necessarily contain every component of the very extended HE absorption complex. The provided rest wavelengths are adequate for coverage checks, not a vetted precision atomic-data list.

The continuum-to-error column is the median `1/error_normalized` within the window, per native 1.3 km/s pixel. It is not the absorption trough's flux/error and not signal-to-noise per resolution element. Resolving power is the nearest 1000-km/s metadata bin's nominal value, not an independently measured line-spread function. Seeing and slit illumination can change the actual line-spread function.

| Target | Fe II rest Å | Nominal observed Å | Valid pixels | Median continuum/error | Nominal R |
|---|---:|---:|---:|---:|---:|
| J051707-441055 | 1608.4509 | 3459.46 | 308 | 127.9 | 62523 |
| J051707-441055 | 2344.2128 | 5041.93 | 307 | 229.7 | 53696 |
| J051707-441055 | 2374.4601 | 5106.99 | 308 | 229.4 | 53696 |
| J051707-441055 | 2382.7642 | 5124.85 | 308 | 238.9 | 53696 |
| J051707-441055 | 2586.6493 | 5563.37 | 308 | 248.4 | 53696 |
| J051707-441055 | 2600.1725 | 5592.45 | 308 | 236.2 | 53696 |
| J034943-381030 | 1608.4509 | 6474.01 | 231 | 86.1 | 53696 |
| J034943-381030 | 2344.2128 | 9435.46 | 230 | 42.6 | 53696 |
| J034943-381030 | 2374.4601 | 9557.20 | 231 | 24.9 | 56936 |
| J034943-381030 | 2382.7642 | 9590.63 | 231 | 21.4 | 56936 |
| J034943-381030 | 2586.6493 | 10411.26 | 231 | 1.3 | 47818 |
| J034943-381030 | 2600.1725 | 10465.69 | 231 | 2.1 | 56936 |

All twelve windows lie inside the wavelength array and have no source-invalid pixels in these intervals. This alone is insufficient for scientific usability. In particular, the high-redshift 2586 and 2600 lines fall near 1.04 μm, where continuum/error is only about 1–2: they are poor precision candidates in this spectrum. The three windows near 9435–9591 Å are in a region requiring serious water-vapour/telluric assessment; sky residuals and blending are also not excluded. The high-redshift 1608 line near 6474 Å is much more promising for initial inspection but still requires contamination checks.

**Inference for design:** the available spectra support immediate analysis-pipeline development and injection/recovery. They do not by themselves provide a clean six-line high-redshift differential measurement. Matching common transitions, adequate negative/positive sensitivity contrast, clean profiles, and exposure-level calibration must precede any cosmological regression. A good low-redshift pilot plus a high-redshift example is not an unbiased redshift sample. A selection function and prospective masking rule are needed before population claims.

## Higher-quality and complementary existing archives

1. **ESPRESSO HE 0515−4414:** the authors publicly release the final spectrum, 17 reduced exposure files, and profile fits. The actual paper demonstrates R≈145000, 16.1 hours and continuum S/N≈105 per 0.4-km/s pixel. Its published alpha constraint is not a Riemann-cutoff constraint. The final FITS is 39,666,240 bytes and was not downloaded in this two-spectrum budget; the entry and exact byte size were verified in the pinned repository tree. [Paper](https://arxiv.org/abs/2112.05819), [data DOI](https://doi.org/10.5281/zenodo.5512490), [pinned final FITS](https://raw.githubusercontent.com/MTMurphy77/ESPRESSO_HE0515-4414/09ee2cb32fb57ace9cda2943e997930bc7813366/Final_spectrum/hes0515m4414.fits). This is an especially useful next baseline and an instrument cross-check on the downloaded SQUAD target.
2. **KODIAQ DR2 / Keck HIRES:** public continuum-normalized 1D flux/error FITS for 300 quasar sightlines, with intermediate products available. The publication reports 831 spectra, while an archive landing page lists 593 coadded products for DR1/DR2; do not silently conflate these product counts. The common supported object count is 300. DR3 is a distinct ESI product and should not be treated as the same HIRES resolution sample. [Official format and retrieval documentation](https://koa.ipac.caltech.edu/applications/KODIAQ/kodiaqDocument.html), [archive landing page](https://koa.ipac.caltech.edu/Datasets/KODIAQ/), [DR2 primary paper](https://arxiv.org/abs/1707.07905). No KODIAQ spectrum was downloaded for this pilot.

## What these data can and cannot test

With a specified transition-response law they can constrain differential line shifts, and some forms of extra broadening or correlated residuals, across absorber redshifts. The first task is to define a response coefficient for each observed transition. A uniform change shared by all transition frequencies is absorbed by the fitted redshift; it is not separately identifiable without an external anchor.

The data contain transition absorption profiles in multiphase, moving gas. They are not a complete ordered atomic/nuclear energy-level sequence, nor direct observations of mathematical Riemann zeros. A generic extra Gaussian width is degenerate with unresolved velocity structure and turbulence unless a predicted dependence on transition/species and adequate controls breaks that degeneracy. Limits on alpha or widths cannot be renamed universal zero-cutoff limits without an explicit physical mapping. No such mapping is established by downloading these spectra.

## Provenance, license, reproducibility

`data/source_manifest.json` records every fetched URL, date, file size and SHA-256, plus pinned source commit `a0cdc8e7b99f2b01a45d919988af9d60e6447d19`. The immutable GitHub source material and the archived tar bytes are retained. Author spectra-host metadata last changed in 2018. The GitHub repository license says CC BY 4.0; the Data Central archive metadata says Creative Commons Attribution Share-Alike without a version. Both original notices are preserved rather than replacing the spectrum-host notice with the repository notice. Cite Murphy et al. (2019) and the data repository in any use.

Run `python code/fetch_public_spectra.py` from this directory (or supply its full path). Cached files are reused; `--refresh` re-downloads. The script never executes downloaded source code and reads the named FITS archive member in memory without extracting arbitrary paths. Dependencies are numpy, requests and astropy. For this session, astropy 8.0.1, pyerfa 2.0.1.5 and astropy-iers-data were installed only in local `.deps/`; the shared Python environment was not upgraded. A regular virtual environment is preferable for redistribution.

Checks completed: catalogue count; target names/redshifts; 9-channel FITS format; monotonic vacuum wavelength grid; endpoint/header agreement; positive-error source-valid mask; full-window geometry; native pixel dispersion; date metadata; per-window error statistics; per-window nominal resolution. No cosmological parameter estimate or detection claim was produced.
