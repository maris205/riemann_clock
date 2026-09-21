# Source-specific rights and attribution

This package preserves the licences and attribution of its source releases.
It does not place all spectra under a new or blanket licence. The Hugging Face
card therefore uses `license: other`; consult the source notices below and the
per-file `source_family` and `license_notice` entries in the manifest/catalogues.

| Material | Source notice | Retained evidence |
| --- | --- | --- |
| ESPRESSO HE 0515−4414 public exposure and final spectrum products, with associated source material | Creative Commons Attribution 4.0 International (CC BY 4.0) | `licenses/ESPRESSO_HE0515-4414-CC-BY-4.0.txt` |
| UVES SQUAD DR1 final-spectrum archives | Creative Commons Attribution Share-Alike; the archived spectrum portal does not specify a version | `provenance/data/raw/catalog_portal_metadata.json` |
| UVES SQUAD DR1 supporting GitHub metadata/code release | Creative Commons Attribution 4.0 International (CC BY 4.0) | `licenses/UVES_SQUAD_DR1-metadata-CC-BY-4.0.txt` |
| Study processed arrays and metadata derived from those spectra | Retain the applicable upstream source obligations; this mirror makes no blanket relicensing grant | Per-file source mapping in `manifest.json` |

The ESPRESSO notice is copied unchanged from
[LICENSE.txt at the pinned source commit](https://github.com/MTMurphy77/ESPRESSO_HE0515-4414/blob/09ee2cb32fb57ace9cda2943e997930bc7813366/LICENSE.txt).
The SQUAD metadata notice is copied unchanged from
[LICENSE at the pinned metadata commit](https://github.com/MTMurphy77/UVES_SQUAD_DR1/blob/a0cdc8e7b99f2b01a45d919988af9d60e6447d19/LICENSE).
The official [CC BY 4.0 licence](https://creativecommons.org/licenses/by/4.0/)
applies to material released under those notices, including its attribution
requirements and warranty disclaimer.

The archived [SQUAD spectrum portal record](https://dmc.datacentral.org.au/api/3/action/package_show?id=uves-squad-dr1)
records `license_id: cc-by-sa`, `license_title: Creative Commons Attribution
Share-Alike`, and `license_url: http://www.opendefinition.org/licenses/cc-by-sa`.
It contains no licence version. This package preserves that uncertainty and
does **not** label those spectra CC BY-SA 4.0 or substitute the metadata/code
repository's CC BY 4.0 notice. No versioned Share-Alike legal text has been
invented for this mirror.

The 52 downloaded spectrum files are redistributed byte for byte. Processed
files are explicitly study derivatives, and their transformations and source
mapping are described in the dataset card and retained project metadata.
Preserve attribution and any applicable Share-Alike obligations when reusing
or distributing those products. Original licence warranty and liability
disclaimers remain in the copied texts; the mirror implies no endorsement by
the original authors or observatories.

Credit and cite the original source authors, not just this mirror:

- Murphy et al. (2022), *A&A* **658**, A123,
  [10.1051/0004-6361/202142257](https://doi.org/10.1051/0004-6361/202142257);
  ESPRESSO data [10.5281/zenodo.5512490](https://doi.org/10.5281/zenodo.5512490).
- Murphy, Kacprzak, Savorgnan & Carswell (2019), *MNRAS* **482**, 3458–3479,
  [10.1093/mnras/sty2834](https://doi.org/10.1093/mnras/sty2834);
  SQUAD data [10.5281/zenodo.1345974](https://doi.org/10.5281/zenodo.1345974).

Study code, context, and the analysis base revision are available at
[maris205/riemann_clock](https://github.com/maris205/riemann_clock/tree/321ae1cd48d69fe4d1091326167a8aaa65f2e597).
