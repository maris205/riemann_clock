#!/usr/bin/env python3
"""Stage the curated spectra dataset without downloading or uploading anything.

Uses only the standard library. Source bytes are never changed. Existing source
checksums are mandatory for every downloaded spectrum; all staged files receive
SHA-256 checksums and project-relative restoration paths in manifest.json.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = Path("experiments/highz_feasibility_2026-09-20")
EXP = ROOT / EXPERIMENT
BASE_COMMIT = "321ae1cd48d69fe4d1091326167a8aaa65f2e597"
REPOSITORY = "https://github.com/maris205/riemann_clock"
DATASET = "dnagpt/riemann-clock-spectra"
EXPECTED_SOURCE = (52, 1_366_293_902)
EXPECTED_PROCESSED = (424, 260_049_012)
ESPRESSO = "ESPRESSO_HE0515-4414"
SQUAD = "UVES_SQUAD_DR1"
ESP_COMMIT = "09ee2cb32fb57ace9cda2943e997930bc7813366"
SOURCE_GLOBS = {
    "data/raw/exposures/*.fits": 17,
    "data/raw/espresso_null/*.fits": 1,
    "data/raw/*.tar.gz": 2,
    "data/raw/archive_expansion/*.tar.gz": 32,
}
PROVENANCE = [
    "data/source_manifest.json",
    "data/raw/README.md",
    "data/raw/Notes_FITS_Files.txt",
    "data/raw/DR1_quasars_master.csv",
    "data/raw/DR1_DLAs.csv",
    "data/raw/catalog_portal_metadata.json",
    "data/raw/espresso_README.md",
    "data/raw/espresso_github_tree.json",
    "data/raw/exposures/manifest.json",
    "data/raw/exposures/provenance_sources.json",
    "data/raw/exposures/provenance_mask_review.json",
    "data/raw/exposures/provenance_upl_line_mask_actions.json",
    "data/raw/espresso_null/manifest.json",
    "data/raw/espresso_null/spectrum_manifest.json",
    "data/raw/espresso_null/followup_sources.json",
    "data/raw/espresso_null/hes0515m4414.upl",
    "data/raw/espresso_null/atmomask.dat",
    "results/archive_expansion/source_manifest.json",
    "results/archive_expansion/strict_source_manifest.json",
    "results/archive_expansion/complete_source_manifest.json",
    "results/archive_expansion/selection_policy.json",
    "results/archive_expansion/complete_scope_policy.json",
    "results/archive_expansion/catalogue_screen.json",
    "results/archive_expansion/strict_catalogue_screen.json",
    "results/noise_covariance/provenance.json",
    "results/noise_covariance/controls.json",
    "results/noise_covariance/intervals.csv",
]
LICENSE_FILES = {
    "data/raw/LICENSE": "licenses/UVES_SQUAD_DR1-metadata-CC-BY-4.0.txt",
    "data/raw/espresso_null/LICENSE.txt": "licenses/ESPRESSO_HE0515-4414-CC-BY-4.0.txt",
}
CSV_FIELDS = [
    "dataset_path", "project_relative_path", "bytes", "sha256", "category",
    "source_family", "target", "instrument", "product_type", "license_notice",
    "source_url", "source_doi", "source_paper_doi", "source_commit",
    "provenance", "source_dataset_paths",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(relative: str):
    return json.loads((EXP / relative).read_text(encoding="utf-8"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def provenance_path(relative: str) -> str:
    return f"provenance/{relative}"


def source_path(relative: str) -> str:
    return "source/" + str(Path(relative).relative_to("data/raw"))


def source_records() -> list[dict]:
    records = {}

    def add(relative, nbytes, checksum, url, manifest, family, target, kind, commit=""):
        row = {
            "dataset_path": source_path(relative),
            "project_relative_path": str(EXPERIMENT / relative),
            "bytes": nbytes, "sha256": checksum, "category": "source_spectrum",
            "source_family": family, "target": target,
            "instrument": "ESPRESSO" if family == ESPRESSO else "UVES",
            "product_type": kind,
            "license_notice": "CC BY 4.0" if family == ESPRESSO else "CC BY-SA; version unspecified by spectrum portal",
            "source_url": url,
            "source_doi": "10.5281/zenodo.5512490" if family == ESPRESSO else "10.5281/zenodo.1345974",
            "source_paper_doi": "10.1051/0004-6361/202142257" if family == ESPRESSO else "10.1093/mnras/sty2834",
            "source_commit": commit,
            "provenance": [provenance_path(manifest)],
        }
        if relative in records:
            old = records[relative]
            require((old["bytes"], old["sha256"], old["source_url"]) == (nbytes, checksum, url),
                    f"Conflicting existing manifests for {relative}")
            old["provenance"].append(provenance_path(manifest))
        else:
            records[relative] = row

    name = "data/raw/exposures/manifest.json"
    for item in read_json(name)["files"]:
        add(item["local_file"], item["bytes"], item["sha256"], item["source_url"],
            name, ESPRESSO, "HE 0515-4414", "extracted_calibrated_exposure", item["commit"])
    name = "data/raw/espresso_null/spectrum_manifest.json"
    item = read_json(name)
    add("data/raw/espresso_null/hes0515m4414.fits", item["bytes"], item["sha256"],
        item["url"], name, ESPRESSO, "HE 0515-4414", "final_coadded_spectrum", ESP_COMMIT)
    name = "data/source_manifest.json"
    for item in read_json(name)["sources"]:
        if item["file"].endswith(".tar.gz"):
            add(item["file"], item["bytes"], item["sha256"], item["url"], name, SQUAD,
                Path(item["file"]).name.split("_Final_")[0], "final_coadded_spectrum_archive")
    for name in [f"results/archive_expansion/{prefix}source_manifest.json"
                 for prefix in ["", "strict_", "complete_"]]:
        for item in read_json(name)["sources"]:
            add(item["archive_file"], item["archive_bytes"], item["archive_sha256"],
                item["archive_url"], name, SQUAD, item["target"], "final_coadded_spectrum_archive")
    actual = set()
    for pattern, expected in SOURCE_GLOBS.items():
        found = sorted(EXP.glob(pattern))
        require(len(found) == expected, f"Expected {expected} files for {pattern}, found {len(found)}")
        actual.update(str(path.relative_to(EXP)) for path in found)
    require(actual == set(records), "Source glob inventory and source manifests differ")
    require((len(records), sum(row["bytes"] for row in records.values())) == EXPECTED_SOURCE,
            "Source inventory totals differ from the fixed release scope")
    for relative, row in records.items():
        path = EXP / relative
        require(path.is_file() and not path.is_symlink(), f"Source is not a regular file: {path}")
        require(path.stat().st_size == row["bytes"] and sha256(path) == row["sha256"],
                f"Source does not match original download manifest: {relative}")
    return sorted(records.values(), key=lambda row: row["dataset_path"])


def original_url_index() -> dict:
    result = {}
    for item in read_json("data/source_manifest.json")["sources"]:
        result[item["file"]] = item
    for item in read_json("data/raw/espresso_null/manifest.json")["files"]:
        relative = item["local_path"].split(str(EXPERIMENT) + "/", 1)[-1]
        result[relative] = item
    for item in read_json("data/raw/espresso_null/followup_sources.json"):
        relative = item["path"].split(str(EXPERIMENT) + "/", 1)[-1]
        result[relative] = item
    return result


def metadata_row(relative: str, destination: str, category: str, family: str,
                 originals: dict) -> dict:
    path = EXP / relative
    row = {
        "dataset_path": destination,
        "project_relative_path": str(EXPERIMENT / relative),
        "bytes": path.stat().st_size, "sha256": sha256(path),
        "category": category, "source_family": family,
        "source_url": f"{REPOSITORY}/blob/{BASE_COMMIT}/{EXPERIMENT / relative}",
        "provenance": [],
    }
    if relative in originals:
        original = originals[relative]
        require(row["sha256"] == original["sha256"] and row["bytes"] == original["bytes"],
                f"Provenance file differs from its original manifest: {relative}")
        row["source_url"] = original["url"]
    return row


def processed_records(sources: list[dict], originals: dict) -> list[dict]:
    rows = []
    source_by_target = {row["target"]: row for row in sources if row["source_family"] == SQUAD}
    exposure_map = {Path(item["array_file"]).name: source_path(item["source_file"])
                    for item in read_json("data/processed/exposures/metadata.json")["exposures"]}
    files = sorted(p for p in (EXP / "data/processed").rglob("*") if p.is_file()
                   and "__pycache__" not in p.parts and p.suffix != ".pyc")
    require((len(files), sum(p.stat().st_size for p in files)) == EXPECTED_PROCESSED,
            "Processed inventory totals differ from the fixed release scope")
    for path in files:
        relative = str(path.relative_to(EXP))
        short = path.relative_to(EXP / "data/processed")
        is_espresso = short.parts[0] == "exposures" or short.name == "espresso_noise_controls.npz"
        family = ESPRESSO if is_espresso else SQUAD
        row = metadata_row(relative, "processed/" + str(short), "processed", family, originals)
        row["license_notice"] = "Derived from CC BY 4.0 release" if is_espresso else "Derived from SQUAD: spectrum portal CC BY-SA (version unspecified); metadata repository CC BY 4.0"
        if short.parts[0] == "exposures":
            row["source_dataset_paths"] = ([exposure_map[short.name]] if short.name in exposure_map
                                            else sorted(exposure_map.values()))
            row["source_dataset_paths"].append(source_path("data/raw/espresso_null/hes0515m4414.fits"))
            row["provenance"] = [provenance_path("data/raw/exposures/manifest.json"),
                                 "processed/exposures/metadata.json"]
        elif short.name == "espresso_noise_controls.npz":
            row["source_dataset_paths"] = [source_path("data/raw/espresso_null/hes0515m4414.fits")]
            row["provenance"] = [provenance_path("results/noise_covariance/provenance.json")]
        else:
            target = short.name.split("_", 1)[0].split(".", 1)[0]
            if target in source_by_target:
                source = source_by_target[target]
                row["target"] = target
                row["source_dataset_paths"] = [source["dataset_path"]]
                row["provenance"] = source["provenance"]
            elif short.name == "feii_coverage_audit.csv":
                row["source_dataset_paths"] = [source_by_target[t]["dataset_path"]
                                                for t in ["J034943-381030", "J051707-441055"]]
                row["provenance"] = [provenance_path("data/source_manifest.json")]
            else:
                raise ValueError(f"Processed source family needs explicit mapping: {relative}")
        rows.append(row)
    return rows


def stage_file(source: Path, destination: Path, row: dict, copy: bool) -> None:
    require(not source.is_symlink(), f"Refusing symlink source {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        require(not destination.is_symlink() and destination.stat().st_size == row["bytes"]
                and sha256(destination) == row["sha256"], f"Existing staged file differs: {destination}")
        return
    if copy:
        shutil.copyfile(source, destination)
    else:
        try:
            os.link(source, destination)
        except OSError:
            shutil.copyfile(source, destination)


def write_generated(output: Path, relative: str, content: str) -> dict:
    path = output / relative
    require(not path.is_symlink(), f"Refusing generated symlink {path}")
    # Atomic replacement avoids writing through a pre-existing hardlink.
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("x", encoding="utf-8", newline="") as handle:
        handle.write(content)
    temporary.replace(path)
    return {"dataset_path": relative, "project_relative_path": None,
            "bytes": path.stat().st_size, "sha256": sha256(path),
            "category": "dataset_documentation", "source_family": "riemann_clock"}


def write_csv(output: Path, filename: str, rows: list[dict]) -> dict:
    import io
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: ";".join(value) if isinstance(value, list) else value
                         for key, value in row.items()})
    return write_generated(output, filename, stream.getvalue())


def verify(output: Path) -> dict:
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    files = manifest["files"]
    listed = [row["dataset_path"] for row in files]
    require(len(listed) == len(set(listed)), "Duplicate dataset paths")
    mapped = [row["project_relative_path"] for row in files if row["project_relative_path"]]
    require(len(mapped) == len(set(mapped)), "Duplicate project restoration paths")
    actual = {str(p.relative_to(output)) for p in output.rglob("*") if p.is_file()}
    require(actual == set(listed) | {"manifest.json"}, "Unlisted files or missing files in staging directory")
    for row in files:
        relative = Path(row["dataset_path"])
        require(not relative.is_absolute() and ".." not in relative.parts, "Unsafe dataset path")
        path = output / relative
        require(not path.is_symlink() and path.stat().st_size == row["bytes"]
                and sha256(path) == row["sha256"], f"Staged checksum failure: {relative}")
        for reference in row.get("source_dataset_paths", []):
            require(reference in listed, f"Missing source reference: {reference}")
        for reference in row.get("provenance", []):
            require(reference in listed, f"Missing provenance reference: {reference}")
    for category, expected in [("source_spectrum", EXPECTED_SOURCE), ("processed", EXPECTED_PROCESSED)]:
        subset = [row for row in files if row["category"] == category]
        require((len(subset), sum(row["bytes"] for row in subset)) == expected,
                f"Incorrect {category} totals")
    with (output / "catalogue.csv").open(newline="", encoding="utf-8") as handle:
        require(len(list(csv.DictReader(handle))) == 52, "Catalogue must contain 52 source product rows")
    with (output / "processed_manifest.csv").open(newline="", encoding="utf-8") as handle:
        require(len(list(csv.DictReader(handle))) == 424, "Processed catalogue must contain 424 rows")
    return {"files_including_manifest": len(files) + 1,
            "bytes_including_manifest": sum((output / p).stat().st_size for p in actual),
            "category_counts": dict(sorted(Counter(row["category"] for row in files).items())),
            "source_manifest_sha256": sha256(output / "manifest.json"),
            "all_staged_sha256_verified": True}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "build/huggingface_dataset")
    parser.add_argument("--copy", action="store_true", help="Copy data instead of using hardlinks")
    parser.add_argument("--verify-only", action="store_true", help="Verify an existing staged package")
    args = parser.parse_args()
    output = args.output.resolve()
    require(output.is_relative_to(ROOT / "build") and output != ROOT / "build",
            "Output must be a subdirectory of the project's ignored build/ directory")
    if not args.verify_only:
        sources = source_records()
        originals = original_url_index()
        processed = processed_records(sources, originals)
        rows = sources + processed
        for relative in PROVENANCE:
            family = (ESPRESSO if any(word in relative for word in ["espresso", "exposures", "noise_covariance"])
                      else SQUAD)
            row = metadata_row(relative, provenance_path(relative), "provenance", family, originals)
            rows.append(row)
        for relative, destination in LICENSE_FILES.items():
            family = ESPRESSO if "espresso" in relative else SQUAD
            rows.append(metadata_row(relative, destination, "license", family, originals))
        output.mkdir(parents=True, exist_ok=True)
        for row in rows:
            stage_file(ROOT / row["project_relative_path"], output / row["dataset_path"], row, args.copy)
        for filename in ["README.md", "LICENSE.md"]:
            template = (ROOT / "data/huggingface" / filename).read_text(encoding="utf-8")
            rows.append(write_generated(output, filename, template))
        rows.append(write_csv(output, "catalogue.csv", sources))
        rows.append(write_csv(output, "processed_manifest.csv", processed))
        manifest = {
            "schema_version": 1, "dataset_id": DATASET, "project_repository": REPOSITORY,
            "analysis_base_commit": BASE_COMMIT,
            "path_mapping": "project_relative_path is relative to the riemann_clock repository root; null denotes dataset-only files",
            "scope": "52 unchanged public extracted/coadded spectrum downloads and all 424 study processed files, with selected provenance and original licenses",
            "manifest_excludes_itself": True,
            "source_spectrum_count": EXPECTED_SOURCE[0], "source_spectrum_bytes": EXPECTED_SOURCE[1],
            "processed_count": EXPECTED_PROCESSED[0], "processed_bytes": EXPECTED_PROCESSED[1],
            "files": sorted(rows, key=lambda row: row["dataset_path"]),
        }
        write_generated(output, "manifest.json", json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(verify(output), indent=2))


if __name__ == "__main__":
    main()
