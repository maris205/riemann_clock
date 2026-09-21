"""Transcribe public experimental summaries; never execute upstream code.

2021 input: publisher supplementary PDF, Tables I and II (pp. 3--4).
2026 input: authors' processed experimental NPZ products, immutable Git SHA
a8b0b38202d6f330297212f8755c8df04242f9d6. These are not raw NMR shots.
Requires pdftotext and numpy; run from any working directory.
"""
from pathlib import Path
import csv
import json
import re
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_parenthetical(value):
    match = re.fullmatch(r"(\d+\.\d+)\((\d+)\)", value)
    if match is None:
        raise ValueError(value)
    digits = len(match[1].split(".")[1])
    return float(match[1]), int(match[2]) * 10.0 ** -digits


def extract_he2021():
    text = subprocess.check_output(
        ["pdftotext", "-layout", str(RAW / "he2021_supplement.pdf"), "-"],
        text=True,
    )
    rows = []
    for line in text.splitlines():
        bits = line.split()
        if not bits or not re.fullmatch(r"\d+", bits[0]):
            continue
        if len(bits) not in (5, 6) or not re.fullmatch(r"\d+\.\d{3}", bits[1]):
            continue
        n = int(bits[0])
        omegas = (5, 8, 12, 16) if n <= 29 else (8, 12, 16)
        assert len(bits[2:]) == len(omegas)
        for omega, value in zip(omegas, bits[2:]):
            measured, sigma = parse_parenthetical(value)
            exact_rounded = float(bits[1])
            rows.append(dict(
                zero_index=n,
                omega_over_J=omega,
                exact_ordinate_published_rounded=exact_rounded,
                measured_ordinate=measured,
                sigma_ordinate_published=sigma,
                parenthetical_as_published=value,
                residual_to_published_rounded=measured-exact_rounded,
                relative_residual_to_published_rounded=(measured-exact_rounded)/exact_rounded,
                sigma_rounded_to_zero=(sigma == 0),
                source_table="I" if n <= 29 else "II",
                source_pdf_page=3 if n <= 29 else 4,
            ))
    assert len(rows) == 269, len(rows)
    assert len({(r["zero_index"], r["omega_over_J"]) for r in rows}) == len(rows)
    assert {r["zero_index"] for r in rows} == set(range(1, 81))
    edge = [r for r in rows if r["zero_index"] == 80]
    assert [(r["omega_over_J"], r["sigma_ordinate_published"]) for r in edge] == [(8, 5.76), (12, 0.05), (16, 0.16)]
    write_csv(OUT / "he2021_supplement_zeros.csv", rows)
    summary = {}
    for omega in (5, 8, 12, 16):
        subset = [r for r in rows if r["omega_over_J"] == omega]
        summaries = []
        for low, high in [(1, 29), (30, 79), (70, 79), (80, 80)]:
            group = [r for r in subset if low <= r["zero_index"] <= high]
            if not group:
                continue
            summaries.append(dict(
                index_range=[low, high], count=len(group),
                median_sigma=float(np.median([r["sigma_ordinate_published"] for r in group])),
                median_absolute_residual=float(np.median([abs(r["residual_to_published_rounded"]) for r in group])),
                max_absolute_residual=float(max(abs(r["residual_to_published_rounded"]) for r in group)),
            ))
        summary[str(omega)] = summaries
    return dict(rows=len(rows), edge_rows=edge, descriptive_summaries=summary,
                warning="Published sigma is an interpolation statistic; not a full error budget. One 189.50(0) entry has rounded-to-zero sigma and cannot receive infinite likelihood weight.")


def extract_wei2026():
    rows = []
    expdir = RAW / "wei2026" / "data" / "exp"
    for name in ["tscan_beta03_exp.npz", "tscan_beta05_exp.npz"]:
        with np.load(expdir / name, allow_pickle=False) as data:
            t = data["t_grid"]
            assert t.shape == (5, 11)
            assert np.array_equal(t.T.reshape(-1), data["t_flat"])
            for i in range(5):
                for j in range(11):
                    row = dict(source_file=name, beta=float(data["beta"]), zero_neighborhood_index=i+1,
                               point_index=j, encoded_dimensionless_time=float(t[i, j]))
                    for key in ["x_mean", "x_std", "y_mean", "y_std", "abs_mean", "abs_std"]:
                        row[key] = float(data[key][i, j])
                    rows.append(row)
    assert len(rows) == 110
    write_csv(OUT / "wei2026_processed_tscan.csv", rows)
    beta_rows = []
    with np.load(expdir / "beta_scan_exp.npz", allow_pickle=False) as data:
        for i, beta in enumerate(data["beta"]):
            row = dict(beta=float(beta), encoded_dimensionless_time=float(data["t_fixed"]))
            for key in ["x_mean", "x_std", "y_mean", "y_std", "abs_mean", "abs_std"]:
                row[key] = float(data[key][i])
            beta_rows.append(row)
    write_csv(OUT / "wei2026_processed_beta_scan.csv", beta_rows)
    reported = [14.12, 20.96, 25.09, 30.44, 32.93]
    write_csv(OUT / "wei2026_published_zero_estimates.csv", [dict(
        zero_index=i+1, reported_ordinate=value,
        uncertainty_not_reported_for_root=True,
        source="Wei et al. (2026), Results: Experimental realization; doi:10.1038/s41467-026-74935-8",
    ) for i, value in enumerate(reported)])
    return dict(tscan_rows=len(rows), beta_scan_rows=len(beta_rows),
                product_level="Authors' processed experimental means and standard deviations; original shots and full cross-point covariance are not present in this downloaded subset.",
                root_estimates="Five root estimates are transcribed from the published Results paragraph; no root uncertainty has been invented.")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    result = dict(he2021=extract_he2021(), wei2026=extract_wei2026())
    (ROOT / "reports" / "experimental_data_extraction.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
