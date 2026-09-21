# Data provenance and meaning

`raw/experimental_sources_download_manifest.json` identifies downloaded files by URL, retrieval date, hash, and where applicable immutable Git commit. Publisher supplements and HTML snapshots are source material, not newly collected observations. The Wei2026 `.npz` files contain authors' processed means and standard deviations; the included subset does not supply all original shots or full covariance matrices.

`processed/he2021_supplement_zeros.csv` preserves all 269 source table entries and their published rounded reference column. `processed/he2021_precision_diagnostics.csv` adds independent mathematical references and neighboring-gap diagnostics. The source reference discrepancies at indices 22 and 23 remain traceable. A zero rounded interpolation uncertainty is flagged, not used as an infinite weight.

`processed/mathematical_reference_zeros.csv` contains 35-digit numerical references at indices 1–81, 4200, and 4201. They are generated with `mpmath.zetazero`; they are not interval-certified results. The validator recomputes them at higher precision.

The Wei2026 time-scan CSV has 110 rows, its inverse-temperature scan has 18 rows, and its published-root CSV has five rows. Those are three different data products; do not add their counts as independent zero measurements. Root uncertainties are not reported. Separate scans are not pooled into a common floor estimate.

`source_project/source_manifest.json` records the unchanged local source PDF and two notebooks, with full archived copies for portable verification. Code-only notebook extracts are archived as `.py.txt` for auditing; the new workflow does not execute them. All new analysis is in `../code`.
