#!/usr/bin/env python3
"""Fetch and audit two public UVES SQUAD DR1 spectra; no spectral fit.

Uses cached files unless --refresh is passed. Read-only HTTPS requests only.
Archives are read in memory; no untrusted archive paths are extracted.
Dependencies: requests, numpy, astropy. Optional local installation: .deps/.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if (ROOT / ".deps").exists():
    sys.path.insert(0, str(ROOT / ".deps"))
import numpy as np
import requests
from astropy.io import fits

RAW, PROCESSED = ROOT / "data/raw", ROOT / "data/processed"
GIT_COMMIT = "a0cdc8e7b99f2b01a45d919988af9d60e6447d19"
BASE = f"https://raw.githubusercontent.com/MTMurphy77/UVES_SQUAD_DR1/{GIT_COMMIT}/"
URLS = {name: BASE + name for name in [
    "README.md", "LICENSE", "Notes_FITS_Files.txt", "DR1_quasars_master.csv"]}
URLS.update({
    "DR1_DLAs.csv": BASE + "DLAs/DR1_DLAs.csv",
    "paper_accepted_2018-10-16.tex": BASE + "Paper/src/paper_accepted_2018-10-16.tex",
    "catalog_portal_metadata.json": "https://dmc.datacentral.org.au/api/3/action/package_show?id=uves-squad-dr1",
    "eso_release_announcement.html": "https://www.eso.org/sci/publications/announcements/sciann17192.html",
    "espresso_README.md": "https://raw.githubusercontent.com/MTMurphy77/ESPRESSO_HE0515-4414/09ee2cb32fb57ace9cda2943e997930bc7813366/README.md",
    "espresso_github_tree.json": "https://api.github.com/repos/MTMurphy77/ESPRESSO_HE0515-4414/git/trees/09ee2cb32fb57ace9cda2943e997930bc7813366?recursive=1",
    "kodiaq_documentation.html": "https://koa.ipac.caltech.edu/applications/KODIAQ/kodiaqDocument.html",
    "J034943-381030_Final_Spectrum.tar.gz": "https://dmc.datacentral.org.au/dataset/c4bad1eb-7e80-49c8-981e-ce44654f79ca/resource/b1ddfba2-62cf-4dbd-9b82-4460fb454688/download/j034943-381030finalspectrum.tar.gz",
    "J051707-441055_Final_Spectrum.tar.gz": "https://dmc.datacentral.org.au/dataset/c4bad1eb-7e80-49c8-981e-ce44654f79ca/resource/c4f30c35-51bc-4db8-acd9-5f9ab5f446c3/download/j051707-441055finalspectrum.tar.gz",
})
TARGETS = {
    "J051707-441055": {"alias": "HE 0515-4414", "z_abs_reference": 1.1508,
                        "z_abs_source": "https://doi.org/10.1093/mnras/stw2543", "half_window_km_s": 200},
    "J034943-381030": {"alias": "Q0347-383 / Q0347-3819", "z_abs_reference": 3.025,
                        "z_abs_source": "SQUAD DR1 master catalogue DLAzabs", "half_window_km_s": 150},
}
FEII_REST_AA = [1608.4509, 2344.2128, 2374.4601, 2382.7642, 2586.6493, 2600.1725]


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def jsonable(value):
    if isinstance(value, bytes):
        return value.decode("ascii", errors="replace").strip()
    if isinstance(value, np.ndarray):
        return [jsonable(x) for x in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def parse_spectrum(name, catalog_row):
    archive = RAW / f"{name}_Final_Spectrum.tar.gz"
    with tarfile.open(archive, "r:gz") as tar:
        member = tar.getmember(f"{name}/{name}.fits")
        if not member.isfile() or member.size > 30_000_000:
            raise ValueError("Unexpected FITS archive member")
        data = tar.extractfile(member).read()
    with fits.open(io.BytesIO(data), memmap=False) as hdul:
        a, h = np.asarray(hdul[0].data, dtype=np.float64), hdul[0].header
        assert a.ndim == 2 and a.shape[0] == 9
        assert int(h["DC-FLAG"]) == 1
        pixel_one_based = np.arange(a.shape[1]) + 1
        wavelength = 10 ** (h["CRVAL1"] + (pixel_one_based - h["CRPIX1"]) * h["CD1_1"])
        assert np.all(np.diff(wavelength) > 0)
        assert np.isclose(wavelength[0], h["UP_WLSRT"], rtol=1e-10)
        assert np.isclose(wavelength[-1], h["UP_WLEND"], rtol=1e-10)
        status = np.rint(a[4]).astype(np.int16)
        valid = (status == 1) & np.isfinite(a[0]) & np.isfinite(a[1]) & (a[1] > 0)
        payload = dict(wavelength_vacuum_heliocentric_AA=wavelength,
                       flux_normalized=a[0], error_normalized=a[1],
                       expected_fluctuation_normalized=a[2], continuum_arbitrary_units=a[3],
                       status=status, contributing_pixels_before_clip=a[5],
                       contributing_pixels_after_clip=a[6], chi2_before_clip=a[7],
                       chi2_after_clip=a[8], valid=valid)
        np.savez_compressed(PROCESSED / f"{name}_squad_dr1.npz", **payload)
        tables = {}
        for i, extension in enumerate(hdul[1:], 1):
            tables[str(i)] = {"name": extension.name, "columns": list(extension.columns.names),
                              "rows": [{col: jsonable(row[col]) for col in extension.columns.names}
                                       for row in extension.data]}
        metadata = {
            "target": name, **TARGETS[name], "catalog": catalog_row,
            "fits_sha256": sha256(data), "n_pixels": a.shape[1], "n_valid_pixels": int(valid.sum()),
            "wavelength_vacuum_heliocentric_AA_min_max": [float(wavelength[0]), float(wavelength[-1])],
            "valid_wavelength_min_max_AA": [float(wavelength[valid][0]), float(wavelength[valid][-1])],
            "dispersion_km_s": h["UP_DISP"], "n_exposures": h["UP_NEXP"], "exposure_seconds": h["UP_TEXP"],
            "observation_date_range_from_exposure_info": [min(row["UTDate"] for row in tables["3"]["rows"]), max(row["UTDate"] for row in tables["3"]["rows"])],
            "exposure_count_note": "Use UP_NEXP for unique observation count; the ExposureInfo table repeats records for detector arms and extracted spectra",
            "hdu_primary_header": dict(h), "tables": tables,
            "source_product": "Continuum-normalized coadded reduced spectrum, not raw photon data",
            "processing": "Only FITS decoding, wavelength reconstruction and validity mask; no resampling, re-normalization or fitted line centroids",
            "covariance": "Neighbouring pixels can be correlated by original redispersion; full pixel covariance is not supplied here",
            "telluric_mask": "No additional telluric, sky, saturation or blend mask applied; valid flag is not proof of uncontaminated absorption",
        }
    (PROCESSED / f"{name}_metadata.json").write_text(json.dumps(metadata, indent=2, default=str) + "\n")
    diagnostics = []
    for rest in FEII_REST_AA:
        centre = rest * (1 + TARGETS[name]["z_abs_reference"])
        velocity = 299792.458 * (wavelength / centre - 1)
        in_window = np.abs(velocity) <= TARGETS[name]["half_window_km_s"]
        good = in_window & valid
        n_window, n_good = int(in_window.sum()), int(good.sum())
        nearest_info = min(tables["1"]["rows"], key=lambda row: abs(row["Wavelength"] - centre))
        diagnostics.append({
            "target": name, "z_abs_reference": TARGETS[name]["z_abs_reference"],
            "transition": f"FeII_{rest:.4f}", "rest_wavelength_AA": rest,
            "nominal_observed_wavelength_AA": centre,
            "half_window_km_s": TARGETS[name]["half_window_km_s"],
            "n_pixels_in_window": n_window, "n_valid_pixels": n_good,
            "masked_fraction_within_array": 1 - n_good/n_window if n_window else None,
            "median_continuum_to_error_per_native_pixel": float(np.median(1/a[1, good])) if n_good else None,
            "median_normalized_flux": float(np.median(a[0, good])) if n_good else None,
            "nearest_metadata_bin_wavelength_AA": nearest_info["Wavelength"],
            "nearest_metadata_nominal_resolving_power": nearest_info["NomResolPower"],
            "nearest_metadata_arc_resolving_power": nearest_info["ArcResolPower"],
            "wavelength_array_contains_full_window": bool(wavelength[0] <= centre * (1 - TARGETS[name]["half_window_km_s"]/299792.458) and wavelength[-1] >= centre * (1 + TARGETS[name]["half_window_km_s"]/299792.458)),
            "interpretation": "Geometric coverage only; identification, blends, tellurics and usable profile need independent screening",
        })
        if n_window:
            np.savez_compressed(PROCESSED / f"{name}_FeII_{rest:.4f}_window.npz",
                                velocity_km_s=velocity[in_window], **{k:v[in_window] for k,v in payload.items()})
    return metadata, diagnostics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Re-download cached sources (read-only remote operations)")
    args = parser.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    manifest_path = ROOT / "data/source_manifest.json"
    previous_manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    previous_entries = {s["file"]: s for s in previous_manifest.get("sources", [])}
    sources = []
    for name, url in URLS.items():
        path = RAW / name
        refreshed = args.refresh or not path.exists()
        if refreshed:
            response = requests.get(url, timeout=(20, 90))
            response.raise_for_status()
            if len(response.content) > 30_000_000:
                raise ValueError("Unexpectedly large source")
            path.write_bytes(response.content)
        content = path.read_bytes()
        previous = previous_entries.get(str(path.relative_to(ROOT)), {})
        digest = sha256(content)
        if not refreshed and previous.get("sha256", digest) != digest:
            raise ValueError(f"Cached source checksum changed: {name}")
        sources.append({"file": str(path.relative_to(ROOT)), "url": url,
                        "bytes": len(content), "sha256": digest,
                        "access_date_utc": datetime.now(timezone.utc).isoformat() if refreshed else
                            previous.get("access_date_utc", "2026-09-20; original pilot retrieval"),
                        "github_commit_if_applicable": GIT_COMMIT if url.startswith(BASE) else None})
    catalog = list(csv.DictReader(io.StringIO((RAW/"DR1_quasars_master.csv").read_text())))
    lookup = {r["Name_Adopt"]: r for r in catalog}
    all_diagnostics, summaries = [], []
    for name in TARGETS:
        meta, diagnostics = parse_spectrum(name, lookup[name])
        all_diagnostics.extend(diagnostics)
        summaries.append({k:v for k,v in meta.items() if k not in ["hdu_primary_header", "tables"]})
    with (PROCESSED / "feii_coverage_audit.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(all_diagnostics[0]))
        writer.writeheader(); writer.writerows(all_diagnostics)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "Data-readiness pilot only; no cosmic-time, alpha-drift or Riemann-cutoff fit",
        "catalog_rows": len(catalog), "final_spectra_count_reported_by_release": 467,
        "catalog_final_spectra_count_by_positive_NumExp": sum(float(r["NumExp"] or 0) > 0 for r in catalog),
        "github_version": GIT_COMMIT,
        "sources": sources, "spectra": summaries,
        "license_note": "GitHub metadata/code repository LICENSE says CC BY 4.0; spectrum-host CKAN portal says CC Attribution Share-Alike (version not specified). Preserve both notices and cite Murphy et al. (2019), DOI 10.1093/mnras/sty2834; do not silently relicense spectra under the repository notice.",
        "downloaded_source_bytes": sum(s["bytes"] for s in sources),
        "limitations": ["Two convenience examples, not a population sample", "Emission redshift differs from absorber redshift", "These are historical coadded observations, not 2026 acquisition data", "No complete atomic/nuclear eigenvalue sequence and no direct Riemann-zero measurement", "Rest wavelengths supplied for coverage checks only, not a vetted precision-constant line list"],
    }
    (ROOT / "data/source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"catalog_rows": len(catalog), "downloaded_source_MB": manifest["downloaded_source_bytes"]/1e6,
                      "spectra": [{"name":m["target"], "pixels":m["n_pixels"], "valid":m["n_valid_pixels"]} for m in summaries]}, indent=2))


if __name__ == "__main__":
    main()
