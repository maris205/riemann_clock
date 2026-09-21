#!/usr/bin/env python3
"""Conditional design forecast, NOT an observed fit or an atomic theory.

The six wavelengths only define a realistic wavelength geometry. The response
coefficients are deliberately synthetic; replacing them by independently derived
atomic coefficients and rerunning is required before testing the Riemann idea.
"""
from pathlib import Path
import csv
import json
import platform
import numpy as np
import scipy
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
YEAR = 31557600.0
MPC_KM = 3.0856775814913673e19
H0 = 67.4
OM = 0.315
TSTAR_SEC = 5.391247e-44
REST_A = np.array([1608.4509, 2344.2128, 2374.4601, 2382.7642, 2586.6493, 2600.1725])
K_RAW = np.array([-1.0, 0.3, 1.0, -0.7, 0.5, -0.2])
K = K_RAW / np.sqrt(np.mean(K_RAW**2))
SEED = 20260920


def age_gyr(z, h0=H0, om=OM):
    """Exact flat matter+Lambda background, radiation neglected, age in Gyr."""
    return (2 * MPC_KM / h0 / (3 * np.sqrt(1 - om)) *
            np.arcsinh(np.sqrt((1 - om) / om) / (1 + np.asarray(z))**1.5) /
            YEAR / 1e9)


def excess(z, p=2, h0=H0, om=OM):
    """h_p(t)-1, evaluated without subtracting nearly equal powers."""
    t0 = age_gyr(0, h0, om)
    l0 = np.log(t0 * 1e9 * YEAR / TSTAR_SEC)
    ell = np.log(age_gyr(z, h0, om) / t0)
    return np.expm1(-p * np.log1p(ell / l0))


def shape(z, p=2, h0=H0, om=OM):
    """Unit response amplitude at z=3, zero at z=0, no predicted amplitude."""
    return excess(z, p, h0, om) / excess(3., p, h0, om)


def projected_response(k=K, degree=1, wave=REST_A):
    x = np.log(wave / np.exp(np.mean(np.log(wave))))
    x /= np.ptp(x)
    nuisance = np.column_stack([x**i for i in range(degree+1)])
    q, _ = np.linalg.qr(nuisance, mode="reduced")
    r = k - q @ (q.T @ k)
    return r, float(r @ r)


def design(n, zlo, zhi, sigma_line, bin_floor, degree=1):
    """Six redshift strata. Floor correlated within each stratum, not /sqrt(N).

    sigma_line: independent line-centroid error assumed equal, m/s.
    bin_floor: residual response-mode systematic std in each redshift stratum,
    m/s; independent between the six strata. No cosmological inference implied.
    """
    assert n % 6 == 0
    z = np.repeat(np.linspace(zlo, zhi, 6), n // 6)
    groups = np.repeat(np.arange(6), n // 6)
    _, k2 = projected_response(degree=degree)
    sigma_abs = sigma_line / np.sqrt(k2)
    cov = np.eye(n) * sigma_abs**2 + bin_floor**2 * (groups[:, None] == groups[None, :])
    factor = cho_factor(cov)
    g = shape(z)
    x = np.column_stack([np.ones(n), g])
    precision_x = cho_solve(factor, x)
    params_cov = np.linalg.inv(x.T @ precision_x)
    estimator = params_cov @ precision_x.T
    sigma_b = np.sqrt(params_cov[1, 1])
    # An arbitrary systematic exactly proportional to g cannot be removed by N.
    align_floor = bin_floor
    sigma_b_aligned = np.sqrt(sigma_b**2 + align_floor**2)
    shape_tests = {}
    for p in (1, 3):
        xp = np.column_stack([np.ones(n), shape(z, p)])
        pp = cho_solve(factor, xp)
        fit = np.linalg.solve(xp.T @ pp, xp.T @ cho_solve(factor, g))
        delta = g - xp @ fit
        unit_distance = np.sqrt(max(float(delta @ cho_solve(factor, delta)), 0.))
        shape_tests[str(p)] = {
            "B_mps_for_expected_delta_chi2_25": float(5 / unit_distance),
            "rms_shape_residual_after_offset_and_amplitude": float(np.sqrt(np.mean(delta**2))),
            "note": "Expected Gaussian separation from one fixed competitor; not a discovery significance or 90%-power result."
        }
    summary = {
        "n_independent_sightlines_assumed": n, "z_range": [zlo, zhi],
        "sigma_line_mps_assumed": sigma_line, "six_bin_correlated_floor_mps_assumed": bin_floor,
        "sigma_single_absorber_mode_mps": float(sigma_abs),
        "sigma_B_at_z3_mps": float(sigma_b),
        "null_95_interval_halfwidth_B_mps": float(norm.ppf(.975) * sigma_b),
        "B_mps_for_5sigma_test_90percent_power": float((5 + norm.ppf(.9)) * sigma_b),
        "sigma_B_with_extra_trend_aligned_floor_mps": float(sigma_b_aligned),
        "signal_shape_discrimination": shape_tests,
        "coverage_all_six_minmax_A": [float(REST_A.min()*(1+zlo)), float(REST_A.max()*(1+zhi))],
    }
    return summary, z, cov, estimator, sigma_b


def run():
    for d in ("results", "figures"):
        (BASE / d).mkdir(exist_ok=True, parents=True)
    scenarios = {
        "UVES_archive_planning": (24, 1.0, 3.0, 200., 50.),
        "ESPRESSO_narrow_redshift_planning": (24, 1.4, 2.0, 50., 10.),
        "ANDES_conditional_2035plus": (60, 1.5, 4.0, 30., 5.),
    }
    rng = np.random.default_rng(SEED)
    output = {
        "status": "SYNTHETIC CONDITIONAL DESIGN; no observed parameter estimates",
        "date": "2026-09-20", "seed": SEED,
        "software": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "background": {"H0_km_s_Mpc": H0, "Omega_m": OM, "flat": True, "radiation": False,
                       "t0_Gyr": float(age_gyr(0)), "tstar_s": TSTAR_SEC},
        "synthetic_response": {"rest_vacuum_wavelength_grid_A_lookup_only": REST_A.tolist(),
                               "K": K.tolist(), "K_origin": "arbitrary toy values; NOT atomic coefficients or alpha q values",
                               "norm_definition": "RMS(K)=1; B is unprojected mode amplitude at z=3 in m/s"},
        "response_retained_squared_norm": {}, "scenarios": {}, "monte_carlo": {},
    }
    for degree in (1, 2):
        r, k2 = projected_response(degree=degree)
        output["response_retained_squared_norm"][str(degree)] = {"norm2": k2,
                    "fraction": k2/float(K @ K), "projected_K": r.tolist()}
    # Exact non-identifiability examples (up to machine precision).
    x = np.log(REST_A)
    for label, response in [("common_shift", np.ones(6)), ("log_wavelength_slope", x)]:
        _, k2 = projected_response(k=response)
        assert k2 < 1e-22
        output.setdefault("degenerate_response_norm2", {})[label] = k2
    for name, args in scenarios.items():
        entry, z, cov, est, sigma_b = design(*args)
        quad, *_ = design(*args, degree=2)
        entry["sigma_B_quadratic_calibration_mps"] = quad["sigma_B_at_z3_mps"]
        entry["background_shape_max_fraction_change_70_030"] = float(np.max(
            np.abs(shape(z, h0=70., om=.30) / shape(z) - 1)))
        output["scenarios"][name] = entry
        checks = []
        nmc = 5000
        for signal_snr in (0., 3., 5+norm.ppf(.9)):
            truth = signal_snr * sigma_b
            mean = truth * shape(z)
            obs = mean + rng.normal(size=(nmc, len(z))) @ np.linalg.cholesky(cov).T
            bhat = obs @ est[1]
            residual = (bhat - truth) / sigma_b
            checks.append({"injected_B_mps": float(truth), "injected_snr": float(signal_snr),
                           "n_simulations": nmc, "mean_bias_in_sigma": float(residual.mean()),
                           "residual_std_in_sigma": float(residual.std(ddof=1)),
                           "coverage_95": float(np.mean(np.abs(residual) < norm.ppf(.975))),
                           "fraction_abs_snr_gt5": float(np.mean(np.abs(bhat/sigma_b) > 5)),
                           "analytic_fraction_abs_snr_gt5": float(norm.sf(5-signal_snr)+norm.cdf(-5-signal_snr))})
        output["monte_carlo"][name] = checks
    ages = []
    for z in (0, 1, 1.1508, 1.4, 2, 3, 3.025, 4, 5):
        ages.append({"z": z, "age_Gyr": float(age_gyr(z)), "h2_minus_one": float(excess(z)),
                     "G2_z3_normalized": float(shape(z))})
    output["time_leverage"] = ages
    with (BASE / "results/forecast_summary.json").open("w") as f:
        json.dump(output, f, indent=2, allow_nan=False)
    with (BASE / "results/time_leverage.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(ages[0])); writer.writeheader(); writer.writerows(ages)
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    zz = np.linspace(0, 5, 300)
    axes[0].plot(zz, 100*excess(zz), color="#1b6e9d", lw=2.5)
    axes[0].set(xlabel="Absorber redshift", ylabel=r"Assumed $h_2(t)-1$ (%)",
                title="Cosmic-time leverage only")
    axes[0].text(.05, .93, r"$h_2=[\ln(t_0/t_*)/\ln(t/t_*)]^2$", transform=axes[0].transAxes, va="top")
    colors = ["#aa5c39", "#428866", "#715ca5"]
    for (name,args), color in zip(scenarios.items(), colors):
        label = name.split("_")[0]
        ns = np.arange(6, 121, 6)
        ss, ideal = [], []
        for n in ns:
            st, *_ = design(int(n), *args[1:])
            nofloor, *_ = design(int(n), *args[1:-1], 0.)
            ss.append(st["sigma_B_at_z3_mps"]); ideal.append(nofloor["sigma_B_at_z3_mps"])
            rows.append({"scenario": name, "N": int(n), "sigma_B_mps": ss[-1], "sigma_B_no_floor_mps": ideal[-1]})
        axes[1].plot(ns, ss, label=label, color=color, lw=2)
        axes[1].plot(ns, ideal, color=color, ls=":", alpha=.7)
    axes[1].set(xlabel="Independent sightlines (assumed)", ylabel=r"Forecast $1\sigma$ on $B$ (m/s)", yscale="log",
                title="Synthetic response / assumed errors")
    axes[1].legend(frameon=False, fontsize=9)
    axes[1].text(.04,.05,"Dotted: no bin-correlated floor",transform=axes[1].transAxes,fontsize=8)
    for p, color in [(1,"#aa5c39"),(3,"#428866")]:
        axes[2].plot(zz[1:], 100*(shape(zz[1:],p)/shape(zz[1:])-1),label=f"p={p} vs p=2",color=color,lw=2)
    axes[2].axhline(0,color="grey",lw=.7)
    axes[2].set(xlabel="Absorber redshift", ylabel="Difference at same z=3 amplitude (%)",
                title="Selecting the exponent is harder")
    axes[2].legend(frameon=False,fontsize=9)
    for ax in axes:
        ax.grid(alpha=.18)
    fig.suptitle("Planning calculation — no cosmic signal measured",fontsize=13)
    for suffix in ("png","pdf"):
        fig.savefig(BASE/f"figures/conditional_sensitivity.{suffix}",dpi=190)
    plt.close(fig)
    with (BASE/"results/sample_size_forecast.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print(json.dumps({"time_leverage": ages, "scenarios": output["scenarios"],
                      "response_projection": output["response_retained_squared_norm"]}, indent=2))


if __name__ == "__main__":
    run()
