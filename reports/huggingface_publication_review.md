# Public spectral dataset publication review

Date: 2026-09-21.

Dataset: [dnagpt/riemann-clock-spectra](https://huggingface.co/datasets/dnagpt/riemann-clock-spectra). Public, ungated, with
immutable release `beb28bc35ed6f0182935582f39b6dfa495d3b9c2`. The pinned revision and manifest
SHA-256 are recorded in `data/huggingface_release.json`.

## Contents and provenance

The package contains 510 files / 1,629,979,693 bytes, excluding the Hub's automatic
`.gitattributes`: 52 downloaded source spectrum products, 424 processed files,
27 provenance files, two source-licence copies, four dataset documents/indexes,
and the manifest itself. All source and processed bytes are unchanged.
The source files reproduce the existing original-download manifests. The data
card distinguishes historical extracted/reduced/coadded products from raw
CCD frames or new observations, and retains source-specific licensing notices.
The SQUAD portal's version-unspecified Share-Alike notice is preserved without
inventing a version or substituting the metadata repository's licence.

## Completed checks

- Independent package review: exact inventory, source hashes/URLs, safe unique
  restoration paths, both CSV indexes, processed-to-source mappings and copied
  licences agree with the retained project records.
- Anonymous remote audit: **514/514 PASS**. All 509 manifest-listed files match:
  440 via the server's SHA-256 metadata and 69 via actual anonymous downloads
  and SHA-256 calculations. The manifest itself matches the pinned SHA-256.
- Fresh anonymous restoration: the 39,666,240-byte ESPRESSO coadd was downloaded
  into an empty temporary destination and exactly matched the original file.
- Local restoration check: all 476 source/processed files, 1,626,342,914 bytes,
  agree with the public manifest. Existing project files were preserved.
- Dataset preview service: both `source_products` and `processed_files` expose
  their `catalogue` split, with no pending or failed entries at the check.
- Independent downloader review covered corrupt or oversized responses,
  preservation of mismatched files, unsafe paths, concurrent installation,
  verify-only behaviour, pinned manifest checks, and transient connection retries.
- Both manuscripts compiled and passed PDF preflight: spectroscopy 20 pages;
  original hypothesis paper 35 pages. HF URI annotations point to this dataset
  on pages 15 and 32 respectively. Scoped visual checks covered spectroscopy
  pages 15–16 and original pages 32–33, with no clipping or overfull boxes.
- Focused manuscript delivery: **72/72 PASS**. Original paper delivery:
  **239/239 PASS**. Scientific source sections, fitting code and numerical
  results were unchanged by this data-publication revision.

The upload initially encountered transfer connection failures; the completed
HTTP transfer and subsequent checks above establish the published result.
The local and remote audits do not rerun the scientific fits or add physical
or external peer-review evidence. Exact checksums and scopes are retained in
`huggingface_release_verification.json` and
`huggingface_download_verification.json`.

## Reproduction

From a current project checkout:

```bash
python code/fetch_huggingface_spectra.py
```

The standard-library command downloads the 52 large source products at the
pinned revision, checks their hashes and restores their original project paths.
Add `--include-processed` for the 424 processed files, or `--verify-only` to
check existing files. Original-source download alternatives remain documented
in `reports/repository_contents.md`.

No authentication is needed to download the public dataset. Upload credentials
were supplied only to the upload process and are absent from the project,
dataset package and Git history; temporary upload/session artifacts are removed
after verification.
