"""Plot the authors' processed 2026 NMR data, without fitting new roots.

Run from any directory with Python, numpy, mpmath and matplotlib installed.
Input CSVs are rebuilt by extract_public_experiments.py. The upstream plotting
scripts were inspected, not imported or executed. abs_mean and abs_std retain
the authors' normalization and reported standard-deviation convention.
"""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"
COMMIT = "a8b0b38202d6f330297212f8755c8df04242f9d6"


def read_csv(name):
    with (DATA / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def column(rows, name):
    return np.asarray([float(row[name]) for row in rows], dtype=float)


def main():
    FIGURES.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    tscan_name = "wei2026_processed_tscan.csv"
    beta_name = "wei2026_processed_beta_scan.csv"
    root_name = "wei2026_published_zero_estimates.csv"
    tscan = read_csv(tscan_name)
    bscan = read_csv(beta_name)
    roots = read_csv(root_name)
    assert len(tscan) == 110 and len(bscan) == 18 and len(roots) == 5
    assert {float(row["beta"]) for row in tscan} == {0.3, 0.5}
    for dataset in (tscan, bscan):
        for name in ("encoded_dimensionless_time", "abs_mean", "abs_std"):
            assert np.isfinite(column(dataset, name)).all(), name
        assert (column(dataset, "abs_std") >= 0).all()

    mp.mp.dps = 45
    gamma_mp = {j: mp.im(mp.zetazero(j)) for j in range(1, 6)}
    gamma = {j: float(value) for j, value in gamma_mp.items()}
    root_rows = []
    for row in roots:
        j = int(row["zero_index"])
        reported = mp.mpf(row["reported_ordinate"])
        root_rows.append({
            "zero_index": j,
            "reference_ordinate_40digits": mp.nstr(gamma_mp[j], 40),
            "reported_ordinate": row["reported_ordinate"],
            "residual_reported_minus_reference": float(reported - gamma_mp[j]),
            "root_standard_deviation": None,
        })

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelsize": 9, "axes.titlesize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.dpi": 220,
    })
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.65))
    fig.subplots_adjust(left=0.083, right=0.985, bottom=0.12,
                        top=0.875, wspace=0.30, hspace=0.49)
    colors = {0.5: "#0072B2", 0.3: "#D55E00"}
    markers = {0.5: "o", 0.3: "s"}
    handles = []
    neighborhood_summary = []
    for j, ax in enumerate(axes.flat[:5], start=1):
        ax.axvline(0, color="0.52", linestyle="--", linewidth=0.8, zorder=1)
        for beta in (0.5, 0.3):
            subset = sorted(
                (row for row in tscan
                 if float(row["beta"]) == beta
                 and int(row["zero_neighborhood_index"]) == j),
                key=lambda row: int(row["point_index"]),
            )
            assert len(subset) == 11
            assert [int(row["point_index"]) for row in subset] == list(range(11))
            t = column(subset, "encoded_dimensionless_time")
            y = column(subset, "abs_mean")
            err = column(subset, "abs_std")
            assert np.allclose(np.diff(t), 0.3, atol=1e-10, rtol=0)
            assert abs(t[5] - gamma[j]) < 1e-6
            h = ax.errorbar(
                t - gamma[j], y, yerr=err, fmt=markers[beta],
                linestyle="none", color=colors[beta], markersize=3.3,
                markerfacecolor="white", markeredgewidth=0.9,
                capsize=2.0, elinewidth=0.8, capthick=0.8, zorder=3,
                label=rf"$\beta={beta:.1f}$",
            )
            if j == 1:
                handles.append(h)
            center = int(np.argmin(abs(t - gamma[j])))
            neighborhood_summary.append({
                "zero_index": j, "beta": beta, "points": len(subset),
                "center_time": float(t[center]),
                "center_minus_reference": float(t[center] - gamma[j]),
                "center_abs_mean": float(y[center]),
                "center_abs_std": float(err[center]),
                "interpretation": "A sampled coherence value, not a fitted root or noise-floor estimate.",
            })
        ax.set_xlim(-1.65, 1.65)
        ax.set_ylim(-0.012, 0.595)
        ax.set_xticks([-1.5, 0, 1.5])
        ax.set_yticks([0, 0.2, 0.4])
        ax.set_title(rf"({chr(96+j)}) $j={j}$, $\gamma_j={gamma[j]:.3f}$", loc="left", pad=5)
        ax.set_xlabel(r"Encoded-time offset $u-\gamma_j$", labelpad=3)
        if j in (1, 4):
            ax.set_ylabel(r"Probe coherence $|\mathcal{L}|$")
        ax.grid(axis="y", color="0.9", linewidth=0.6)

    ax = axes.flat[5]
    bscan = sorted(bscan, key=lambda row: float(row["beta"]))
    assert np.ptp(column(bscan, "encoded_dimensionless_time")) == 0
    ax.axvline(0.5, color="0.52", linestyle="--", linewidth=0.8, zorder=1)
    ax.errorbar(column(bscan, "beta"), column(bscan, "abs_mean"),
                yerr=column(bscan, "abs_std"), fmt="o", linestyle="none",
                color="#4C566A", markerfacecolor="white", markersize=3.3,
                capsize=2, elinewidth=0.8, capthick=0.8, zorder=3)
    ax.set(xlim=(0, 0.95), ylim=(-0.005, 0.15), xlabel=r"$\beta$",
           ylabel=r"$|\mathcal{L}|$")
    ax.set_xticks([0, 0.25, 0.5, 0.75])
    ax.set_yticks([0, 0.05, 0.10, 0.15])
    ax.set_title(r"(f) Separate scan at $u\simeq\gamma_1$", loc="left", pad=5)
    ax.grid(axis="y", color="0.9", linewidth=0.6)
    fig.legend(handles, [r"$\beta=0.5$", r"$\beta=0.3$"], loc="upper center",
               bbox_to_anchor=(0.50, 1.005), ncol=2, frameon=False,
               handlelength=1.5, columnspacing=2.2)
    fig.text(0.5, 0.925, "Wei et al. (2026): processed NMR measurements; bars show reported SD",
             ha="center", va="center", fontsize=9)
    fig.savefig(FIGURES / "wei2026_experiment.pdf", metadata={
        "Title": "Wei et al. 2026 processed NMR data near the first five Riemann zeros",
        "Author": "riemann_clock reproducibility analysis",
        "Subject": "Processed experimental means and reported SD; no new root fitting",
        "CreationDate": None, "ModDate": None,
    })
    fig.savefig(FIGURES / "wei2026_experiment.png")
    plt.close(fig)

    with (RESULTS / "wei2026_published_root_residuals.csv").open(
            "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(root_rows[0]))
        writer.writeheader()
        writer.writerows(root_rows)

    table = [r"\begin{tabular}{rrrr}", r"\toprule",
             r"$j$ & Reference $\gamma_j$ & Published estimate & Estimate $-\gamma_j$ \\",
             r"\midrule"]
    for row in root_rows:
        table.append(f"{row['zero_index']} & {float(row['reference_ordinate_40digits']):.9f} & "
                     f"{float(row['reported_ordinate']):.2f} & "
                     f"{row['residual_reported_minus_reference']:+.6f} " + r"\\")
    table.extend([r"\bottomrule", r"\end{tabular}"])
    (RESULTS / "wei2026_root_table.tex").write_text("\n".join(table) + "\n", encoding="utf-8")

    summary = {
        "article_doi": "10.1038/s41467-026-74935-8",
        "publication_date": "2026-07-01", "arxiv_first_public_date": "2025-11-14",
        "repository_commit": COMMIT,
        "data_level": "Authors' processed NMR means and SD, not raw repeated signals or shots.",
        "time_scan_points": len(tscan), "beta_scan_points": len(bscan),
        "plotted_quantity": "abs_mean directly from authors' processed products; no rescaling",
        "error_bars": "abs_std; one reported standard deviation, not standard errors or root uncertainties",
        "time_definition": "u = epsilon_0 * laboratory_time / hbar; not cosmic age",
        "physical_experiment_scope": "Five-qubit NMR; first five zero neighborhoods only",
        "high_zero_scope": "The trillionth-zero results in the source are numerical simulations",
        "reference_calculation": {"method": "mpmath.zetazero", "decimal_precision": mp.mp.dps,
                                  "mpmath_version": mp.__version__},
        "source_csv_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in (tscan_name, beta_name, root_name)
        },
        "root_estimates": root_rows,
        "root_residual_descriptive_summary": {
            "mean_absolute_residual": float(np.mean([abs(row["residual_reported_minus_reference"]) for row in root_rows])),
            "maximum_absolute_residual": float(max(abs(row["residual_reported_minus_reference"]) for row in root_rows)),
            "fit_performed_here": False,
            "no_chi_squared": "Root standard deviations and full covariance are not published in this subset.",
        },
        "sampled_neighborhood_centers": neighborhood_summary,
        "warnings": [
            "The source's rounded polynomial-fit root values are transcribed, not fitted again here.",
            "Means and standard deviations are retained; repeated signals and full covariance are unavailable.",
            "Nonzero sampled minima are not estimates of a universal residual floor.",
            "The beta scan and time scan are separate data series, not repetitions to pool.",
            "Formal publication lies in the last six months; first public preprint does not.",
        ],
    }
    (RESULTS / "wei2026_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"figures": ["wei2026_experiment.pdf", "wei2026_experiment.png"],
                      "time_scan_points": len(tscan), "beta_scan_points": len(bscan),
                      "root_estimates": root_rows}, indent=2))


if __name__ == "__main__":
    main()
