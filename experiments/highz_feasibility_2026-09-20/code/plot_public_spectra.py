#!/usr/bin/env python3
"""Plot the twelve real SQUAD DR1 Fe II windows; no fit or injection.

Requires the processed files made by fetch_public_spectra.py. The source flux,
statistical error and native pixel positions are displayed without smoothing.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed"
OUT = ROOT / "figures"
TARGETS = ["J051707-441055", "J034943-381030"]
REST = [1608.4509, 2344.2128, 2374.4601, 2382.7642, 2586.6493, 2600.1725]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with (DATA / "feii_coverage_audit.csv").open() as stream:
        rows = {(row["target"], float(row["rest_wavelength_AA"])): row
                for row in csv.DictReader(stream)}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": "#273345", "text.color": "#172439",
                         "axes.titleweight": "medium", "pdf.fonttype": 42})
    fig, axes = plt.subplots(6, 2, figsize=(12.2, 14.6), sharex="col")
    colours = ["#126782", "#A14B14"]
    for j, target in enumerate(TARGETS):
        for i, rest in enumerate(REST):
            ax = axes[i, j]
            audit = rows[(target, rest)]
            with np.load(DATA / f"{target}_FeII_{rest:.4f}_window.npz", allow_pickle=False) as d:
                velocity = d["velocity_km_s"]
                good = d["valid"]
                flux = np.where(good, d["flux_normalized"], np.nan)
                error = np.where(good, d["error_normalized"], np.nan)
            ax.fill_between(velocity, flux-error, flux+error, step="mid", color=colours[j], alpha=0.18, linewidth=0)
            ax.step(velocity, flux, where="mid", color=colours[j], linewidth=0.8)
            ax.axhline(1, color="#9CA3AF", linestyle=":", linewidth=0.7)
            ax.axhline(0, color="#D4D8DD", linewidth=0.5)
            ax.axvline(0, color="#7A5195", linestyle="--", linewidth=0.8)
            ax.grid(axis="x", alpha=0.14)
            half = float(audit["half_window_km_s"])
            ax.set_xlim(-half, half)
            low, high = float(np.nanmin(flux-error)), float(np.nanmax(flux+error))
            if j == 1 and i >= 4:
                # Preserve all observed central values and error bands in noisy windows.
                pad = 0.06*(high-low)
                ax.set_ylim(low-pad, high+pad)
                ax.text(.98, .88, "LOW C/N · expanded y-axis", ha="right", va="top",
                        transform=ax.transAxes, fontsize=8, color="#A14B14",
                        bbox={"facecolor":"white", "edgecolor":"none", "alpha":.8, "pad":2})
            else:
                ax.set_ylim(min(-.12, low-.03), max(1.17, high+.06))
            if j == 1 and i in (1, 2, 3):
                ax.text(.98, .12, "Telluric assessment required", ha="right", va="bottom",
                        transform=ax.transAxes, fontsize=8, color="#805032",
                        bbox={"facecolor":"white", "edgecolor":"none", "alpha":.8, "pad":1.5})
            cnr = float(audit["median_continuum_to_error_per_native_pixel"])
            centre = float(audit["nominal_observed_wavelength_AA"])
            ax.set_title(f"Fe II {rest:.1f} Å  →  {centre:.1f} Å   |   C/N {cnr:.1f}",
                         loc="left", fontsize=9, pad=6)
            if j == 0:
                ax.set_ylabel("Normalized flux")
            if i == 5:
                ax.set_xlabel("Velocity relative to nominal absorber redshift (km/s)")
    fig.suptitle("Real archival spectra: twelve nominal Fe II windows", x=.5, y=.981,
                 fontsize=17, fontweight="bold")
    fig.text(.5, .959, "UVES SQUAD DR1 · native pixels and published statistical errors · no model fit or injected signal",
             ha="center", fontsize=10, color="#455468")
    fig.text(.277, .928, "HE 0515−4414   |   z_abs = 1.1508", ha="center", fontsize=12, color=colours[0], fontweight="bold")
    fig.text(.757, .928, "Q0347−3819   |   z_abs = 3.025", ha="center", fontsize=12, color=colours[1], fontweight="bold")
    fig.legend(handles=[Line2D([0],[0], color="#33465E", linewidth=1, label="Observed normalized flux"),
                        Patch(facecolor="#33465E", alpha=.18, label="Published flux ± 1σ statistical error"),
                        Line2D([0],[0], color="#7A5195", linestyle="--", linewidth=.9, label="Nominal absorber redshift")],
               loc="lower center", bbox_to_anchor=(.5,.055), ncol=3, frameon=False, fontsize=9)
    fig.text(.5,.037, "C/N = median(1 / normalized error), per 1.3 km/s pixel; it is not a fitted line-position precision.",
             ha="center", fontsize=9, color="#455468")
    fig.text(.5,.022, "No telluric/blend correction is applied here. The high-z windows near 1.04 μm have C/N ≈ 1–2 and are poor precision candidates.",
             ha="center", fontsize=9, color="#455468")
    fig.subplots_adjust(left=.075, right=.975, top=.903, bottom=.112, hspace=.39, wspace=.18)
    fig.savefig(OUT / "public_feii_windows.png", dpi=180, facecolor="white")
    fig.savefig(OUT / "public_feii_windows.pdf", facecolor="white")
    plt.close(fig)
    print("Saved figures/public_feii_windows.png and .pdf (real observed flux, no fit).")


if __name__ == "__main__":
    main()
