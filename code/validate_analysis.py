"""Independent numerical and transcription checks for the descriptive reanalysis.

Reads existing products without importing or running either analysis module.
Uses Decimal parsing and 50-digit zero evaluations as cross-checks. Numerical
zero evaluations are not certified interval proofs. Writes an internal check
report, not a peer-review or author-approval attestation.
"""
from collections import Counter
from decimal import Decimal
from pathlib import Path
import csv
import hashlib
import json
import math
import re
import statistics
import subprocess

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 50


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def close(actual, expected, *, atol=1e-11, rtol=1e-12):
    assert math.isclose(float(actual), float(expected), abs_tol=atol, rel_tol=rtol), (actual, expected)


def main():
    checks = {}
    # Hash each independently downloaded input and the unchanged original files.
    manifests = [
        (ROOT / "data/source_project/source_manifest.json", None),
        (ROOT / "data/raw/experimental_sources_download_manifest.json", ROOT / "data/raw"),
    ]
    matched = 0
    original_files_checked = 0
    for path, base in manifests:
        for row in json.loads(path.read_text()):
            original = ROOT / "data/source_project" / row["archived_path"] if base is None else base / row["path"]
            contents = original.read_bytes()
            assert hashlib.sha256(contents).hexdigest() == row["sha256"], original
            assert len(contents) == row["bytes"], original
            matched += 1
            if base is None and Path(row["source"]).exists():
                assert hashlib.sha256(Path(row["source"]).read_bytes()).hexdigest() == row["sha256"]
                original_files_checked += 1
    checks["source_file_hashes_matched"] = matched
    checks["original_external_files_unchanged_if_available"] = original_files_checked

    # Re-parse the two table pages directly, avoiding the extraction code/parser.
    pdftext = subprocess.check_output(
        ["pdftotext", "-layout", str(ROOT / "data/raw/he2021_supplement.pdf"), "-"], text=True)
    pattern = re.compile(r"^\s*(\d+)\s+(\d+\.\d{3})\s+((?:\d+\.\d+\(\d+\)\s*){3,4})$", re.M)
    source = {}
    for index_text, exact_text, measured_text in pattern.findall(pdftext):
        index = int(index_text)
        values = re.findall(r"(\d+\.\d+)\((\d+)\)", measured_text)
        omegas = [5, 8, 12, 16] if len(values) == 4 else [8, 12, 16]
        for omega, (center_text, uncertainty_text) in zip(omegas, values):
            center = Decimal(center_text)
            sigma = Decimal(uncertainty_text).scaleb(center.as_tuple().exponent)
            source[index, omega] = (Decimal(exact_text), center, sigma)
    assert len(source) == 269
    rows = read_csv(ROOT / "data/processed/he2021_supplement_zeros.csv")
    assert len(rows) == 269
    zero_sigma_keys = []
    for row in rows:
        key = int(row["zero_index"]), int(row["omega_over_J"])
        exact, center, sigma = source[key]
        assert Decimal(row["exact_ordinate_published_rounded"]) == exact
        assert Decimal(row["measured_ordinate"]) == center
        close(row["sigma_ordinate_published"], sigma, atol=1e-14)
        if sigma == 0:
            zero_sigma_keys.append(key)
            assert row["sigma_rounded_to_zero"] == "True"
        else:
            assert row["sigma_rounded_to_zero"] == "False"
    assert zero_sigma_keys == [(74, 8)]
    assert [source[80, omega][2] for omega in (8, 12, 16)] == [Decimal("5.76"), Decimal("0.05"), Decimal("0.16")]
    checks["he2021_transcription"] = {"matched_records": len(rows), "records_by_drive": dict(Counter(int(r["omega_over_J"]) for r in rows)),
        "rounded_zero_sigma_keys": zero_sigma_keys, "index80_sigmas": [5.76, .05, .16]}

    # Regenerate references at higher precision; independently check small zeta residuals.
    refs = read_csv(ROOT / "data/processed/mathematical_reference_zeros.csv")
    gamma = {}
    max_difference = mp.mpf(0)
    max_zeta_residual = mp.mpf(0)
    for row in refs:
        index = int(row["zero_index"])
        computed = mp.im(mp.zetazero(index))
        gamma[index] = computed
        difference = abs(computed - mp.mpf(row["ordinate_decimal"]))
        max_difference = max(max_difference, difference)
        assert difference < mp.mpf("1e-30")
        zeta_residual = abs(mp.zeta(mp.mpc(mp.mpf(".5"), computed)))
        assert zeta_residual < mp.mpf("1e-43")
        max_zeta_residual = max(max_zeta_residual, zeta_residual)
    assert gamma[19] < 76 < gamma[20]
    assert mp.im(mp.zetazero(3681)) < 4200 < mp.im(mp.zetazero(3682))
    checks["mathematical_references"] = {"entries_recomputed": len(gamma), "precision_digits": 50,
        "max_absolute_difference": str(max_difference), "max_zeta_residual": str(max_zeta_residual),
        "gamma80": str(gamma[80]), "gamma81": str(gamma[81]), "N76": 19, "N4200": 3681,
        "certified_interval_proof": False}

    diagnostics = read_csv(ROOT / "data/processed/he2021_precision_diagnostics.csv")
    reconstructed = {}
    bad_source_refs = set()
    for row in diagnostics:
        j, omega = int(row["zero_index"]), int(row["omega_over_J"])
        exact, center, sigma = source[j, omega]
        neighbors = [gamma[j + 1] - gamma[j]]
        if j > 1:
            neighbors.append(gamma[j] - gamma[j - 1])
        gap = min(neighbors)
        residual = mp.mpf(str(center)) - gamma[j]
        for field, value in [("reference_ordinate", gamma[j]), ("reference_neighbor_gap", gap),
                             ("residual", residual), ("absolute_residual_over_neighbor_gap", abs(residual) / gap),
                             ("published_sigma_over_neighbor_gap", mp.mpf(str(sigma)) / gap)]:
            close(row[field], value)
        reconstructed[j, omega] = (residual, mp.mpf(str(sigma)), gap)
        if abs(mp.mpf(str(exact)) - gamma[j]) > mp.mpf(".0005"):
            bad_source_refs.add(j)
    assert bad_source_refs == {22, 23}
    group_rows = read_csv(ROOT / "results/he2021_group_summaries.csv")
    for row in group_rows:
        omega, lo, hi = (int(row[k]) for k in ("omega_over_J", "first_index", "last_index"))
        selected = [values for (j, om), values in reconstructed.items() if om == omega and lo <= j <= hi]
        assert len(selected) == int(row["count"])
        close(row["median_absolute_residual"], statistics.median(abs(v[0]) for v in selected))
        close(row["median_published_sigma"], statistics.median(v[1] for v in selected))
        assert int(row["central_residual_exceeds_one_gap"]) == sum(abs(r) > gap for r, s, gap in selected)
        assert int(row["published_sigma_exceeds_one_gap"]) == sum(s > gap for r, s, gap in selected)
    beyond18 = {omega: [j for (j, om), (r, s, gap) in reconstructed.items() if om == omega and j > 18 and abs(r) > gap]
                for omega in (8, 12, 16)}
    assert {omega: len(js) for omega, js in beyond18.items()} == {8: 3, 12: 1, 16: 2}
    checks["he2021_diagnostics"] = {"matched_rows": len(diagnostics), "matched_group_summaries": len(group_rows),
        "source_reference_discrepancy_indices": sorted(bad_source_refs), "beyond18_residual_exceeds_gap_indices": beyond18,
        "index80_neighbor_gap": str(min(gamma[80] - gamma[79], gamma[81] - gamma[80]))}

    scans = read_csv(ROOT / "data/processed/wei2026_processed_tscan.csv")
    seen = set()
    numeric_fields = ["x_mean", "x_std", "y_mean", "y_std", "abs_mean", "abs_std"]
    for filename in {row["source_file"] for row in scans}:
        with np.load(ROOT / "data/raw/wei2026/data/exp" / filename, allow_pickle=False) as raw:
            for row in (r for r in scans if r["source_file"] == filename):
                i, j = int(row["zero_neighborhood_index"]) - 1, int(row["point_index"])
                key = filename, i, j
                assert key not in seen
                seen.add(key)
                assert float(row["beta"]) == float(raw["beta"])
                assert float(row["encoded_dimensionless_time"]) == raw["t_grid"][i, j]
                for field in numeric_fields:
                    assert float(row[field]) == raw[field][i, j]
    betas = read_csv(ROOT / "data/processed/wei2026_processed_beta_scan.csv")
    with np.load(ROOT / "data/raw/wei2026/data/exp/beta_scan_exp.npz", allow_pickle=False) as raw:
        for i, row in enumerate(betas):
            assert float(row["beta"]) == raw["beta"][i]
            assert float(row["encoded_dimensionless_time"]) == float(raw["t_fixed"])
            for field in numeric_fields:
                assert float(row[field]) == raw[field][i]
    assert len(scans) == 110 and len(betas) == 18
    roots = read_csv(ROOT / "data/processed/wei2026_published_zero_estimates.csv")
    assert len(roots) == 5
    assert all(row["uncertainty_not_reported_for_root"] == "True" for row in roots)
    checks["wei2026_data_integrity"] = {"time_scan_rows": len(scans), "beta_scan_rows": len(betas),
        "numeric_cells_exactly_matched_to_npz": 8 * (len(scans) + len(betas)),
        "root_uncertainties_explicitly_missing": True, "upstream_code_executed": False}

    # Check inverse-spacing identity, resource monotonicity, and floor saturation.
    # This checks the conditional illustration, not any real instrument bound.
    scenarios = read_csv(ROOT / "results/illustrative_resource_scenarios.csv")
    by_floor = {}
    for row in scenarios:
        m, b, delta, logheight, q = (mp.mpf(row[k]) for k in
            ("shots", "assumed_floor", "assumed_delta", "mean_density_envelope_log10_T", "q"))
        close(delta**2 - b**2, 1 / m)
        logheight_over_2pi = logheight * mp.log(10) - mp.log(2 * mp.pi)
        mean_spacing = 2 * mp.pi / logheight_over_2pi
        close(delta, q * mean_spacing)
        assert q == mp.mpf(".25")
        by_floor.setdefault(b, []).append((m, logheight))
        if b > 0:
            assert logheight < mp.log10(2 * mp.pi) + 2 * mp.pi * q / (b * mp.log(10))
    for sequence in by_floor.values():
        assert all(y2 > y1 and m2 > m1 for (m1, y1), (m2, y2) in zip(sequence, sequence[1:]))
    checks["conditional_resource_scenarios"] = {"rows": len(scenarios), "q": .25,
        "identities_and_monotonicity_pass": True, "not_empirical_limits": True}

    summary = json.loads((ROOT / "results/analysis_summary.json").read_text())
    aging = summary["aging_scenario"]
    t0, seconds_per_year, tstar = (mp.mpf(str(aging[k])) for k in ("t0_years", "seconds_per_year", "tstar_seconds"))
    log0 = mp.log(t0 * seconds_per_year / tstar)
    normalized = lambda t, p=2: (log0 / mp.log(t * seconds_per_year / tstar)) ** p
    derivative = mp.diff(normalized, t0)
    halfyear_change = normalized(t0 + mp.mpf(".5")) - 1
    close(aging["d_log_delta_cos_per_year"], derivative, atol=1e-27)
    close(aging["halfyear_fractional_change_delta_cos"], halfyear_change, atol=1e-27)
    age_rows = read_csv(ROOT / "results/conditional_aging_scenarios.csv")
    for row in age_rows:
        close(row["normalized_delta_cos"], normalized(mp.mpf(row["cosmic_age_Gyr"]) * 10**9, int(row["exponent"])))
    checks["conditional_aging"] = {"rows": len(age_rows), "derivative_per_year": str(derivative),
        "halfyear_fractional_change": str(halfyear_change), "age_behavior_not_fitted": True}
    checks["status"] = "PASS: internal numerical and transcription consistency checks"
    output = ROOT / "results/numerical_validation.json"
    output.write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
