"""Reproducible descriptive reanalysis and explicitly conditional illustrations.

No cosmological parameters are fitted. Published interpolation uncertainties are
not treated as complete Gaussian errors. Numerical zero values are high-precision
mpmath references, not certified interval enclosures.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
FIG = ROOT / "figures"
RESULT = ROOT / "results"
PAPER = ROOT / "paper"
mp.mp.dps = 35
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42,
                     "savefig.dpi": 180})
COLORS = {5: "#999999", 8: "#0072B2", 12: "#D55E00", 16: "#009E73"}


def references():
    path = DATA / "mathematical_reference_zeros.csv"
    indices = list(range(1, 82)) + [4200, 4201]
    if path.exists():
        rows = pd.read_csv(path, dtype={"ordinate_decimal": str})
        assert rows.zero_index.tolist() == indices
        assert (rows.precision_decimal_digits == 35).all()
    else:
        records = [{"zero_index": j,
                    "ordinate_decimal": mp.nstr(mp.im(mp.zetazero(j)), 35),
                    "precision_decimal_digits": 35,
                    "method": "mpmath.zetazero", "mpmath_version": mp.__version__}
                   for j in indices]
        rows = pd.DataFrame(records)
        rows.to_csv(path, index=False)
    return {int(r.zero_index): mp.mpf(r.ordinate_decimal) for r in rows.itertuples()}


def main():
    for path in (FIG, RESULT, PAPER):
        path.mkdir(parents=True, exist_ok=True)
    gamma = references()
    source = DATA / "he2021_supplement_zeros.csv"
    df = pd.read_csv(source)
    assert len(df) == 269
    df["reference_ordinate"] = df.zero_index.map(lambda j: float(gamma[j]))
    df["reference_neighbor_gap"] = df.zero_index.map(
        lambda j: float(gamma[2]-gamma[1]) if j == 1 else
        float(min(gamma[j]-gamma[j-1], gamma[j+1]-gamma[j])))
    df["residual"] = df.measured_ordinate - df.reference_ordinate
    df["absolute_residual_over_neighbor_gap"] = abs(df.residual)/df.reference_neighbor_gap
    df["published_sigma_over_neighbor_gap"] = df.sigma_ordinate_published/df.reference_neighbor_gap
    # Source reference column has small transcription/rounding discrepancies.
    # Preserve it verbatim; diagnose against independently computed values.
    reference_discrepancies = df[abs(df.reference_ordinate - df.exact_ordinate_published_rounded) > .0005]
    reference_discrepancies = reference_discrepancies.drop_duplicates("zero_index")
    df.to_csv(DATA / "he2021_precision_diagnostics.csv", index=False)

    stats = []
    for omega, group in df.groupby("omega_over_J", sort=True):
        for lo, hi in ((1, 18), (19, 80), (61, 80), (1, 80)):
            part = group[group.zero_index.between(lo, hi)]
            if part.empty:
                continue
            stats.append(dict(omega_over_J=int(omega), first_index=lo, last_index=hi,
                count=len(part), median_absolute_residual=float(abs(part.residual).median()),
                median_published_sigma=float(part.sigma_ordinate_published.median()),
                central_residual_exceeds_one_gap=int((part.absolute_residual_over_neighbor_gap > 1).sum()),
                published_sigma_exceeds_one_gap=int((part.published_sigma_over_neighbor_gap > 1).sum())))
    pd.DataFrame(stats).to_csv(RESULT / "he2021_group_summaries.csv", index=False)
    rows = [r for r in stats if r["omega_over_J"] != 5 and r["first_index"] in (1,61)
            and not (r["first_index"] == 1 and r["last_index"] == 80)]
    tex = [r"\begin{tabular}{crrrr}", r"\toprule",
           r"$\Omega/J$ & Index range & Count & Median $|e_j|$ & Median $s_j$ \\", r"\midrule"]
    for r in rows:
        tex.append(f"{r['omega_over_J']} & {r['first_index']}--{r['last_index']} & {r['count']} & "
                   f"{r['median_absolute_residual']:.4f} & {r['median_published_sigma']:.3f} " + r"\\")
    tex += [r"\bottomrule", r"\end{tabular}"]
    (PAPER / "he2021_summary_table.tex").write_text("\n".join(tex)+"\n")

    fig, axes = plt.subplots(4, 2, figsize=(10.0, 9.3), sharex=True,
                             gridspec_kw={"width_ratios": [1.35, 1]})
    for row, omega in enumerate((5, 8, 12, 16)):
        g = df[df.omega_over_J == omega]
        ax = axes[row, 0]
        ax.axhline(0, color="0.3", lw=.8)
        ax.errorbar(g.zero_index, g.residual, yerr=g.sigma_ordinate_published,
                    fmt="o", ms=2.7, color=COLORS[omega], elinewidth=.65, capsize=1.4)
        ax.set_ylabel(r"$\hat\gamma_j-\gamma_j$")
        ax.text(.02, .93, rf"$\Omega/J={omega}$", transform=ax.transAxes, va="top")
        ax.margins(y=.18)
        # All error bars remain visible; axis ranges are intentionally independent.
        ax = axes[row, 1]
        ax.plot(g.zero_index, g.absolute_residual_over_neighbor_gap, "o", ms=2.8,
                color=COLORS[omega], label=r"$|e_j|/d_j$")
        positive = g.published_sigma_over_neighbor_gap > 0
        ax.plot(g[positive].zero_index, g[positive].published_sigma_over_neighbor_gap,
                "x", ms=3.5, color="0.25", label=r"$s_j/d_j$")
        ax.axhline(1, color="0.4", ls="--", lw=.8)
        ax.set_yscale("log")
        ax.set_ylabel("Ratio to local gap")
        ax.set_ylim(1e-4, 20)
        if row == 0:
            ax.legend(loc="lower right", frameon=False, fontsize=8)
        if omega == 8:
            ax.annotate(r"$s_{74}=0$ (rounded; omitted)", xy=(.03,.07),
                        xycoords="axes fraction", fontsize=7.5)
    for ax in axes[-1]:
        ax.set_xlabel("Zero index $j$")
        ax.set_xlim(0, 82)
    axes[0, 0].set_title("Central residuals and published interpolation errors")
    axes[0, 1].set_title("Descriptive neighboring-gap comparison")
    fig.tight_layout()
    for suffix in ("pdf", "png"):
        fig.savefig(FIG / f"he2021_precision.{suffix}", bbox_inches="tight")
    plt.close(fig)

    year_seconds = 365.25*86400
    t0_years = 13.8e9
    tstar_seconds = 5.391247e-44
    log0 = float(np.log(t0_years*year_seconds/tstar_seconds))
    drift = -2/(t0_years*log0)
    # log1p/expm1 avoid cancellation for sub-year age changes.
    halfyear_fractional = float(np.expm1(-2*np.log1p(np.log1p(.5/t0_years)/log0)))
    checks = {
        "method": "Descriptive analysis; no cosmic fit, no joint Gaussian likelihood",
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "rows": len(df), "unique_zero_indices": 80,
        "published_reference_discrepancies_over_half_last_digit": reference_discrepancies[
            ["zero_index", "exact_ordinate_published_rounded", "reference_ordinate"]].to_dict(orient="records"),
        "gamma_18": str(gamma[18]), "gamma_19": str(gamma[19]),
        "gamma_80": str(gamma[80]), "gamma_4200": str(gamma[4200]), "gamma_4201": str(gamma[4201]),
        "N_76": int(mp.nzeros(76)), "N_4200": int(mp.nzeros(4200)),
        "old_log_formula_lab": float(np.log(1e10)/(np.pi*np.sqrt(.01))),
        "old_log_formula_cosmic": float(np.log(1e60)/(np.pi*np.sqrt(.0001))),
        "old_quadratic_assumption_wait_ordinate_4200_to_4201_years": t0_years*np.expm1(2*np.log1p(1/4200)),
        "old_quadratic_assumption_wait_index_4200_to_4201_years": float(t0_years*((gamma[4201]/gamma[4200])**2-1)),
        "wait_values_are_not_physical_forecasts": True,
        "aging_scenario": {"t0_years": t0_years, "tstar_seconds": tstar_seconds,
            "seconds_per_year": year_seconds, "log_t0_over_tstar": log0,
            "p": 2, "d_log_delta_cos_per_year": drift,
            "halfyear_fractional_change_delta_cos": halfyear_fractional,
            "delta_cos0_fitted": False},
        "edge_estimates": df[df.zero_index==80].to_dict(orient="records"),
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
            "pandas": pd.__version__, "mpmath": mp.__version__, "matplotlib": matplotlib.__version__},
    }
    (RESULT / "analysis_summary.json").write_text(json.dumps(checks, indent=2)+"\n")

    # All values below are conditional scenarios, not fitted experimental results.
    m = np.logspace(0, 4, 161)
    q = .25
    scenario_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.05))
    for floor, color in ((0, "#0072B2"), (.1, "#D55E00"), (.3, "#009E73")):
        delta = np.sqrt(1/m + floor**2)  # kappa*tau*visibility = 1; illustrative attainment
        log10_height = np.log10(2*np.pi) + 2*np.pi*q/(delta*np.log(10))
        axes[0].plot(m, log10_height, color=color, label=rf"Assumed floor $b={floor:g}$")
        scenario_rows.extend(dict(shots=float(mm), assumed_floor=floor,
            assumed_delta=float(dd), mean_density_envelope_log10_T=float(tt),
            q=q, kappa_tau_visibility=1) for mm, dd, tt in zip(m,delta,log10_height))
    pd.DataFrame(scenario_rows).to_csv(RESULT / "illustrative_resource_scenarios.csv", index=False)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Independent repetitions $M$")
    axes[0].set_ylabel(r"Illustrative $\log_{10} T_{\rm env}$")
    axes[0].set_title("(a) Mean-density envelope under assumed precision")
    axes[0].legend(frameon=False, fontsize=8)
    age_gyr = np.geomspace(1, 100, 201)
    aging_rows = []
    for p, color, ls in ((1,"#999999","--"),(2,"#0072B2","-"),(3,"#D55E00",":")):
        normalized = (log0/np.log(age_gyr*1e9*year_seconds/tstar_seconds))**p
        axes[1].plot(age_gyr, normalized, color=color, ls=ls, label=f"Exponent $p={p}$")
        aging_rows.extend(dict(cosmic_age_Gyr=float(t), exponent=p,
            normalized_delta_cos=float(y)) for t,y in zip(age_gyr, normalized))
    pd.DataFrame(aging_rows).to_csv(RESULT / "conditional_aging_scenarios.csv", index=False)
    axes[1].axvline(13.8, color="0.5", lw=.8, ls="--")
    axes[1].axhline(1, color="0.7", lw=.7)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Assumed cosmic age (Gyr)")
    axes[1].set_ylabel(r"$\delta_{\rm c}(t)/\delta_{\rm c}(t_0)$")
    axes[1].set_title("(b) Conditional aging; normalization only")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    for suffix in ("pdf", "png"):
        fig.savefig(FIG / f"conditional_resources.{suffix}", bbox_inches="tight")
    plt.close(fig)
    print(json.dumps({"rows": len(df), "gamma_80": str(gamma[80]),
                      "aging_fractional_drift_per_year": drift,
                      "outputs": "processed diagnostics, summaries, 2 figures"}, indent=2))


if __name__ == "__main__":
    main()
