# Repository contents and recovery of large source files

The Git repository includes both manuscript PDFs and LaTeX sources, figures,
analysis code, processed spectra, saved fit arrays and numerical endpoints,
validation reports, source manifests, and upstream attribution/license files.
Scientific results include unsuccessful or capped optimizer endpoints where
these are part of the recorded analysis history.

Build directories, local Python dependencies and caches are ignored. The large
third-party downloads below are also ignored; their source URLs and checksums
remain in the supplied manifests. The original local files are preserved.

## Restore from the public Hugging Face snapshot

The [public dataset](https://huggingface.co/datasets/dnagpt/riemann-clock-spectra)
contains all 52 source spectrum products listed below, plus 424 processed files.
These are historical extracted/reduced/coadded products and their derivatives;
they are not new observations or original detector frames.

From the project root, run:

```bash
python code/fetch_huggingface_spectra.py
```

This standard-library command requires no token. It reads the immutable dataset
revision and manifest hash from [the release record](../data/huggingface_release.json),
checks each file's size and SHA-256, and restores files to their original project
paths. Existing matching files are retained; a mismatch stops without overwriting.
Unlike re-running the original fetch/processing scripts, this restore command does
not regenerate source manifests or processed arrays. Use `--include-processed` to
also restore the 424 processed files, `--verify-only` to check existing files, or
`--dry-run` to list the required paths and total size.

## Original-source recovery alternatives

| Files omitted from Git | Recovery script in the spectroscopy experiment |
| --- | --- |
| 17 ESPRESSO S2D exposure FITS files | `code/exposure_fetch.py` |
| ESPRESSO combined-spectrum FITS file | `code/fetch_espresso_null_sources.py` |
| Two initial SQUAD spectrum archives | `code/fetch_public_spectra.py` |
| 32 additional SQUAD spectrum archives | `code/archive_expand_all.py` |

These omissions do not prevent compiling either paper from its saved figures
and tables. Full raw-data validation and analyses that reopen original FITS or
archives require the relevant downloads. Install the experiment's Python
dependencies and run the following commands from the repository root, choosing
only the source families needed for the intended analysis:

```bash
cd experiments/highz_feasibility_2026-09-20
python code/fetch_espresso_null_sources.py
python code/exposure_fetch.py
python code/fetch_public_spectra.py
python code/archive_expand_all.py
```

The ESPRESSO scripts verify the preserved SHA-256 values or pinned Git blob
identities. The SQUAD scripts use the preserved source URLs and catalogue
selection; compare newly downloaded bytes against the supplied source manifests.
Retrieval can update access timestamps, manifests and processed arrays. Preserve
the supplied snapshots (or use a separate checkout) before running those scripts;
review `git diff` afterwards. A changed hash needs investigation and must not be
silently accepted as the original input.

The retained sources include third-party material under its original terms;
this repository does not assign a new blanket license to those materials.
The laboratory source subset under `data/raw/` remains included and is separate
from the larger astronomy downloads listed above.
