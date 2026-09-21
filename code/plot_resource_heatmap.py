"""Heatmap of the paper's conditional resource model; no empirical anchors.

The color is log10 of an ordinate envelope, not a measured cutoff or zero count.
All parameters are illustrative. Neither grid values nor colors are clipped.
"""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
Q = 0.25


def log_height(repetitions, floor):
    delta = np.hypot(1 / np.sqrt(repetitions), floor)
    return np.log10(2 * np.pi) + 2 * np.pi * Q / (delta * np.log(10))


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "pdf.fonttype": 42, "savefig.dpi": 240})
    for sub in ["figures", "results"]:
        (ROOT / sub).mkdir(exist_ok=True)
    repetitions = np.geomspace(1, 1e4, 321)
    floors = np.geomspace(.03, 1, 241)
    mesh_m, mesh_b = np.meshgrid(repetitions, floors)
    heights = log_height(mesh_m, mesh_b)
    assert np.isfinite(heights).all()
    assert (np.diff(heights, axis=1) > 0).all()
    assert (np.diff(heights, axis=0) < 0).all()
    np.savez_compressed(ROOT / "results/resource_heatmap_grid.npz",
        repetitions=repetitions, assumed_floor=floors, log10_ordinate_envelope=heights,
        q=Q, kappa_tau_visibility=1.)

    fig, ax = plt.subplots(figsize=(10.8, 6.7))
    fig.subplots_adjust(left=.1, right=.88, bottom=.14, top=.86)
    im = ax.pcolormesh(mesh_m, mesh_b, heights, cmap="viridis", shading="gouraud",
                       vmin=float(heights.min()), vmax=float(heights.max()), rasterized=True)
    contours = ax.contour(mesh_m, mesh_b, heights,
        levels=[2, 3, 5, 10, 15, 20], colors="white", linewidths=.7, alpha=.75)
    label_m = 6000.
    label_positions = [(label_m, float(np.sqrt(
        (2*np.pi*Q/((level-np.log10(2*np.pi))*np.log(10)))**2 - 1/label_m)))
        for level in [2, 3, 5, 10, 15, 20]]
    ax.clabel(contours, inline=True, fmt="%g", fontsize=8, manual=label_positions)
    # A precision floor dominates where b > M^(-1/2); equality is a crossover.
    crossover_m = np.geomspace(1.01, 1 / .03**2, 200)
    ax.plot(crossover_m, 1 / np.sqrt(crossover_m), ls="--", color="#ffdf70", lw=1.5)
    ax.text(13, .195, r"$b=M^{-1/2}$", color="#ffdf70", fontsize=10,
            path_effects=[pe.withStroke(linewidth=2.5, foreground="#25334b")])
    samples = []
    for name, m, b, offset, align in [
        ("A", 100., .3, (-14, 15), "right"),
        ("B", 3000., .3, (12, 12), "left"),
        ("C", 3000., .05, (12, 10), "left"),
    ]:
        value = float(log_height(m, b))
        samples.append({"label": name, "repetitions": m, "assumed_floor": b,
                        "log10_ordinate_envelope": value, "ordinate_envelope": float(10**value)})
        ax.scatter(m, b, marker="*", s=220, facecolor="#ffcf56", edgecolor="white",
                   linewidth=1.1, zorder=5)
        ax.annotate(name, (m,b), xytext=offset, textcoords="offset points", color="white",
                    fontsize=13, fontweight="bold", ha=align,
                    path_effects=[pe.withStroke(linewidth=2, foreground="#29354e")])
    ax.annotate("", xy=(2000,.3), xytext=(170,.3),
                arrowprops={"arrowstyle":"->", "color":"white", "lw":1.3})
    ax.annotate("", xy=(3000,.065), xytext=(3000,.235),
                arrowprops={"arrowstyle":"->", "color":"white", "lw":1.3})
    ax.text(.20, .89, "Precision floor dominates", transform=ax.transAxes,
            color="white", fontsize=12, fontweight="medium")
    ax.text(.055, .085, "Statistical term dominates", transform=ax.transAxes,
            color="white", fontsize=11)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(repetitions[0], repetitions[-1])
    ax.set_ylim(floors[0], floors[-1])
    ax.set_xlabel("Independent repetitions $M$", labelpad=8)
    ax.set_ylabel("Assumed precision floor $b$ (ordinate units)", labelpad=8)
    ax.grid(which="major", color="white", ls=":", alpha=.14)
    cbar = fig.colorbar(im, ax=ax, pad=.026, fraction=.045)
    cbar.set_label(r"$\log_{10} T_{\rm env}$", labelpad=9)
    cbar.set_ticks([2,5,10,15,20])
    fig.suptitle("A resource landscape for Riemann-zero estimation", y=.967,
                 fontsize=17, fontweight="semibold")
    fig.text(.5,.905, r"Conditional mean-density model: $\delta^2=M^{-1}+b^2$, $q=1/4$",
             ha="center", fontsize=11, color="#45516a")
    fig.text(.1,.035, "Stars A–C: illustrative scenarios  •  Dashed line: equal error contributions",
             fontsize=9.5, color="#45516a")
    for suffix in ["pdf", "png"]:
        fig.savefig(ROOT / f"figures/resource_precision_heatmap.{suffix}", bbox_inches="tight")
    plt.close(fig)
    summary = {"interpretation":"Conditional mean-density envelope; not a measured cutoff or zero count",
        "shape":list(heights.shape), "q":Q, "kappa_tau_visibility":1,
        "repetitions_range":[float(repetitions[0]),float(repetitions[-1])],
        "floor_range":[float(floors[0]),float(floors[-1])],
        "log10_ordinate_envelope_range":[float(heights.min()),float(heights.max())],
        "values_clipped":False, "empirical_anchors":False, "samples":samples}
    (ROOT / "results/resource_heatmap_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
