---
pretty_name: Riemann Clock Public Spectra and Processed Products
language:
- en
license: other
license_name: source-specific-licenses
license_link: LICENSE.md
tags:
- astronomy
- spectroscopy
- quasar
- reproducibility
size_categories:
- n<1K
configs:
- config_name: source_products
  default: true
  data_files:
  - split: catalogue
    path: catalogue.csv
- config_name: processed_files
  data_files:
  - split: catalogue
    path: processed_manifest.csv
---

# Riemann Clock spectra

This dataset accompanies the exploratory spectroscopy study in
[maris205/riemann_clock](https://github.com/maris205/riemann_clock), with analysis
base commit [`321ae1cd48d69fe4d1091326167a8aaa65f2e597`](https://github.com/maris205/riemann_clock/tree/321ae1cd48d69fe4d1091326167a8aaa65f2e597).
It mirrors public, previously reduced, extracted, or coadded quasar spectra and
the study's processed arrays. It contains no newly acquired observations and no
raw CCD frames. The project directory name `data/raw/` means downloaded inputs.
This exploratory study reports no detection of Riemann-clock physics, varying
fundamental constants, or a cosmological time signal.

| Contents | Files | Bytes |
| --- | ---: | ---: |
| ESPRESSO extracted, wavelength-calibrated exposure FITS | 17 | 1,108,454,400 |
| ESPRESSO final combined FITS | 1 | 39,666,240 |
| UVES SQUAD individual final-spectrum archives | 34 | 218,173,262 |
| Study processed arrays, metadata, and coverage table | 424 | 260,049,012 |

The 52 unchanged source products total 1,366,293,902 bytes; the source and
processed products together total 1,626,342,914 bytes. Provenance, licence
copies, catalogues, and this card add a small amount of storage. Counts refer
to files, not independent observations or independent measurements: the
ESPRESSO combined spectrum and its 17 exposures reuse the same photons, and
many processed windows reuse pixels from the same coadd. The UVES subset
contains 34 sightlines, including HE 0515−4414 also observed with ESPRESSO.

## Download and restore

From a current [project checkout](https://github.com/maris205/riemann_clock),
restore all 52 source products to their original project locations with:

```bash
python code/fetch_huggingface_spectra.py
```

That recovery script uses the dataset revision and manifest checksum pinned in
`data/huggingface_release.json`, verifies file hashes, and preserves the original
source manifests. Add `--include-processed` to restore the processed files too,
or use `--verify-only` to check an existing local restoration.

Install `huggingface_hub`, then download the complete dataset:

```python
from huggingface_hub import snapshot_download

snapshot_download(
    repo_id="dnagpt/riemann-clock-spectra",
    repo_type="dataset",
    local_dir="riemann-clock-spectra",
)
```

For repeatable downloads, supply `revision="<dataset commit SHA>"` after selecting
a published dataset revision. The Hugging Face viewer previews the 52-row
`catalogue.csv`; the `processed_files` configuration previews a separate 424-row
index. FITS, NPZ, and tar archives remain downloadable binary files. The CSV
rows are file metadata, not one row per wavelength pixel.

`manifest.json` lists each distributed file except itself, with `dataset_path`,
`project_relative_path`, byte count, SHA-256, category, and source family. Source
products also carry their original URL, citation DOI, and existing provenance
manifest references. The two CSV indexes carry the same source and processed
file hashes. `source_dataset_paths` in processed records connects each product
to its source download; several files have more than one source.

The dataset uses short `source/` and `processed/` paths. To restore the study's
layout, copy each file with a non-null `project_relative_path` to that path
**relative to the cloned `riemann_clock` repository root**, checking its size
and SHA-256 first. For example,
`source/espresso_null/hes0515m4414.fits` maps to
`experiments/highz_feasibility_2026-09-20/data/raw/espresso_null/hes0515m4414.fits`.
Keep existing files when their hashes agree and reject mismatches rather than
silently replacing local changes. Dataset-only cards and indexes have null
restoration paths. Additional study code, fit inputs, results, and documentation
are available in the linked GitHub repository; this is a curated spectrum
package, not a replacement for the entire project checkout.

## Sources and attribution

**ESPRESSO HE 0515−4414:** the 17 extracted exposures and final combined spectrum
come from the [authors' repository](https://github.com/MTMurphy77/ESPRESSO_HE0515-4414/tree/09ee2cb32fb57ace9cda2943e997930bc7813366),
pinned at `09ee2cb32fb57ace9cda2943e997930bc7813366`; data DOI
[10.5281/zenodo.5512490](https://doi.org/10.5281/zenodo.5512490).
Cite Murphy et al. (2022), “Fundamental physics with ESPRESSO: Precise limit on
variations in the fine-structure constant towards the bright quasar
HE 0515−4414,” *Astronomy & Astrophysics* **658**, A123,
[10.1051/0004-6361/202142257](https://doi.org/10.1051/0004-6361/202142257),
[arXiv:2112.05819](https://arxiv.org/abs/2112.05819).
The source release uses CC BY 4.0. Its original licence and README are retained.

**UVES SQUAD DR1:** the 34 final-spectrum `.tar.gz` downloads come from the
[Data Central spectrum portal](https://dmc.datacentral.org.au/dataset/uves-squad-dr1),
whose archived metadata retains its earlier Swinburne hosting URLs. The exact
download URL for each file is in `catalogue.csv`. The supporting
[metadata/code repository](https://github.com/MTMurphy77/UVES_SQUAD_DR1/tree/a0cdc8e7b99f2b01a45d919988af9d60e6447d19)
is pinned at `a0cdc8e7b99f2b01a45d919988af9d60e6447d19`; data DOI
[10.5281/zenodo.1345974](https://doi.org/10.5281/zenodo.1345974).
Cite Murphy, Kacprzak, Savorgnan & Carswell (2019), “The UVES Spectral Quasar
Absorption Database (SQUAD) Data Release 1: The first 10 million seconds,”
*Monthly Notices of the Royal Astronomical Society* **482**(3), 3458–3479,
[10.1093/mnras/sty2834](https://doi.org/10.1093/mnras/sty2834),
[arXiv:1810.06136](https://arxiv.org/abs/1810.06136).
The spectrum portal specifies Creative Commons Attribution Share-Alike
**without a version number**. The metadata/code repository separately provides
CC BY 4.0. These are distinct notices; this mirror does not relicense the
spectra under the metadata repository's licence. See [LICENSE.md](LICENSE.md).

## Processed products and limitations

The 424 processed files are preserved exactly as present in the study:

- `processed/archive_expansion/`: 32 UVES coadd NPZ files, 32 metadata JSON files,
  and 324 absorption-window NPZ files for 36 absorber entries.
- The root `processed/J034943-381030_*` and `processed/J051707-441055_*` files:
  two UVES spectra, their metadata, and 12 Fe II windows. The latter target is
  HE 0515−4414, but these files come from UVES, not ESPRESSO.
- `processed/feii_coverage_audit.csv`: the initial two-target UVES coverage audit.
- `processed/exposures/`: 17 ESPRESSO native-pixel segment NPZ files and exposure
  metadata, with the study's masking and count scaling recorded there.
- `processed/espresso_noise_controls.npz`: ESPRESSO coadd control windows for
  the study's noise diagnostics.

The producing scripts are retained in the project's
[experiment code directory](https://github.com/maris205/riemann_clock/tree/321ae1cd48d69fe4d1091326167a8aaa65f2e597/experiments/highz_feasibility_2026-09-20/code):
`fetch_public_spectra.py`, `archive_fetch.py`, `archive_expand_all.py`,
`exposure_analysis.py`, and `empirical_noise_controls.py`.

NPZ files use NumPy named arrays; inspect their array names and the accompanying
metadata before analysis. Original FITS headers and upstream FITS format notes
are retained. Source and processed wavelength conventions differ by product
(for example, native ESPRESSO vacuum barycentric wavelengths and published
UVES vacuum heliocentric coadds); follow each product's metadata. Processed
masks reflect the specific study choices and do not guarantee an uncontaminated
or precision-ready line sample. Coaddition, interpolation, blends, continuum
choices, and shared photons can induce correlations. Complete pixel covariance
and raw detector calibration material are not supplied by this package.

The selected archives are a study subset of SQUAD DR1, not the entire release
or a statistically complete population sample. The archived source manifests
and catalogues can mention upstream files beyond this package; only files in
`manifest.json` are distributed here. No publisher PDFs, general software
manuals, or unrelated analysis scripts are included.

## Rebuild and verify

With the original project data restored, run the repository's
`python code/prepare_huggingface_dataset.py`. It stages this release under the
ignored `build/huggingface_dataset/`, verifies all 52 source hashes against
existing download manifests, and checks every staged file. Re-running the
builder with unchanged inputs produces the same metadata and manifest. Use
`--copy` to stage independent copies; the default uses hardlinks where possible,
so treat staged binary files as immutable. `--verify-only` checks an existing
package. The script does not access the network or upload data.
