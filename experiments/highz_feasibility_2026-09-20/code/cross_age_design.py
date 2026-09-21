#!/usr/bin/env python3
"""Fixed six-target geometry, not a measurement of cosmic atomic variation.

No measured shifts enter these rank, curvature or hypothetical-error forecasts.
Per-target unconstrained differential wavelength/gas terms destroy physical
identifiability. All finite precision values below condition on removing them.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform

import numpy as np
import scipy
from scipy.linalg import solve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/cross_age"
SOURCE = ROOT / "results/archive_expansion/quality_summary.json"
ATOMIC = ROOT / "data/raw/espresso_null/MM_VPFIT_2013-11-10.dat"
YEAR = 31557600.0
MPC_KM = 3.0856775814913673e19
C_M_S = 299792458.0
H0 = 67.4
OM = 0.315
TSTAR = 5.391247e-44


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def age_gyr(z):
    z = np.asarray(z, dtype=float)
    return (2 * MPC_KM / H0 / (3 * np.sqrt(1 - OM)) *
            np.arcsinh(np.sqrt((1 - OM) / OM) / (1 + z)**1.5) /
            (1e9 * YEAR))


def excess(z, p=2, tstar=TSTAR):
    """Stable evaluation with dimensional t/tstar always in seconds."""
    t0 = age_gyr(0)
    l0 = np.log(t0 * 1e9 * YEAR / tstar)
    ell = np.log(age_gyr(z) / t0)
    return np.expm1(-p * np.log1p(ell / l0))


def curve(z, key):
    if key == "inverse_log_p2_tstar_1second":
        return excess(z, 2, 1.) / excess(3., 2, 1.)
    if key.startswith("inverse_log_p"):
        p = int(key[-1])
        return excess(z, p) / excess(3., p)
    if key == "log_cosmic_age":
        return np.log(age_gyr(0) / age_gyr(z)) / np.log(age_gyr(0) / age_gyr(3.))
    if key == "linear_cosmic_age":
        return (age_gyr(0) - age_gyr(z)) / (age_gyr(0) - age_gyr(3.))
    if key == "linear_redshift":
        return np.asarray(z, dtype=float) / 3.
    raise ValueError(key)


def reference_wavelength(key):
    rows = []
    for line in ATOMIC.read_text().splitlines():
        x = line.split()
        if len(x) > 5 and x[0] == "FeII" and key <= float(x[1]) < key + 1:
            rows.append([float(x[1]), float(x[2])])
    a = np.asarray(rows)
    return float(np.sum(a[:, 0] * a[:, 1]) / a[:, 1].sum())


def project(vector, nuisance):
    """Unweighted design, normalized columns; SVD rank is relative 1e-12."""
    nuisance = np.asarray(nuisance, dtype=float)
    scale = np.linalg.norm(nuisance, axis=0)
    normalized = nuisance[:, scale > 0] / scale[scale > 0]
    u, singular, _ = np.linalg.svd(normalized, full_matrices=False)
    rank = int(np.count_nonzero(singular > singular[0] * 1e-12))
    q = u[:, :rank]
    residual = vector - q @ (q.T @ vector)
    return residual, rank, singular


def generate():
    OUT.mkdir(parents=True, exist_ok=True)
    targets = sorted(json.loads(SOURCE.read_text())["readiness_candidates"], key=lambda x: x["z_abs"])
    z = np.array([x["z_abs"] for x in targets])
    g = curve(z, "inverse_log_p2")
    one = np.ones(len(z))
    wav = np.array([reference_wavelength(k) for k in [2374, 2382]])
    separation = np.diff(wav)[0] * (1 + z)
    baseline = g - g.mean()
    full_information = float(baseline @ baseline)
    designs = {
        "common_offset": one[:, None],
        "common_offset_and_shared_wavelength_slope": np.column_stack([one, separation]),
        "one_wavelength_slope_per_target": np.column_stack([one, np.diag(separation)]),
        "one_environmental_pair_offset_per_target": np.column_stack([one, np.eye(len(z))]),
    }
    rank_results = {}
    for name, nuisance in designs.items():
        residual, rank, sv = project(g, nuisance)
        _, augmented_rank, augmented_sv = project(g, np.column_stack([nuisance, g]))
        norm2 = float(residual @ residual)
        structural_zero = augmented_rank == rank
        rank_results[name] = dict(
            nuisance_rank=rank, augmented_rank=augmented_rank,
            identifiable_additional_rank=augmented_rank - rank,
            squared_projected_response_raw=norm2,
            squared_projected_response=0. if structural_zero else norm2,
            information_fraction_vs_common_offset=0. if structural_zero else norm2 / full_information,
            singular_values=sv.tolist(), augmented_singular_values=augmented_sv.tolist(),
            projected_response=residual.tolist(), structural_nonidentifiability=structural_zero)
    keys = ["inverse_log_p1", "inverse_log_p3", "log_cosmic_age", "linear_cosmic_age", "linear_redshift"]
    comparison = {}
    for key in keys:
        alternative = curve(z, key)
        x = np.column_stack([one, alternative])
        beta = np.linalg.lstsq(x, g, rcond=None)[0]
        residual = g - x @ beta
        ac = alternative - alternative.mean()
        rho = float(baseline @ ac / np.linalg.norm(baseline) / np.linalg.norm(ac))
        norm = float(np.linalg.norm(residual))
        comparison[key] = dict(
            centered_correlation=rho,
            fitted_offset=float(beta[0]), fitted_amplitude=float(beta[1]),
            residual_after_offset_and_amplitude=residual.tolist(),
            residual_norm=norm,
            fractional_shape_residual_relative_to_centered_p2=float(norm / np.linalg.norm(baseline)),
            max_abs_residual_for_hypothetical_A_100_m_s=float(np.max(np.abs(residual)) * 100),
            rms_residual_for_hypothetical_A_100_m_s=float(norm / np.sqrt(len(z)) * 100),
            hypothetical_sigma_each_for_expected_delta_chi2_9_at_A_100_m_s=100 * norm / 3,
            warning="Fixed Gaussian mean separation only. A=100 m/s is arbitrary, not fitted or predicted. No detection significance, power calculation, or multiplicity correction.")
    forecasts = []
    for sigma in [10., 30., 100.]:
        for lab_sigma in [0., 20., 100.]:
            cov = np.eye(len(z)) * sigma**2 + lab_sigma**2 * np.ones((len(z), len(z)))
            x = np.column_stack([one, g])
            precision_x = solve(cov, x, assume_a="pos")
            pcov = np.linalg.inv(x.T @ precision_x)
            forecasts.append(dict(sigma_independent_pair_m_s=sigma,
                sigma_shared_laboratory_pair_m_s=lab_sigma,
                conditional_sigma_A_m_s=float(np.sqrt(pcov[1, 1])),
                conditional_sigma_sample_endpoint_difference_m_s=float(np.sqrt(pcov[1, 1]) * np.ptp(g))))
    rows = []
    for i, target in enumerate(targets):
        rows.append(dict(target=target["target"], z_abs=float(z[i]), cosmic_age_Gyr=float(age_gyr(z[i])),
            G2_z3_normalized=float(g[i]), h2_minus_1=float(excess(z[i])),
            observed_2374_AA=float(wav[0] * (1 + z[i])),
            observed_2382_AA=float(wav[1] * (1 + z[i])),
            observed_pair_separation_AA=float(separation[i]),
            velocity_slope_m_s_per_AA_that_exactly_mimics_A_100_m_s=float(100 * g[i] / separation[i])))
    shared_slope_fit = np.linalg.lstsq(np.column_stack([one, separation]), 100 * g, rcond=None)[0]
    shared_slope_residual = 100 * g - np.column_stack([one, separation]) @ shared_slope_fit
    slope_scenarios = []
    for slope_per_1000_AA in [50., 200., 500.]:
        pair_effect = slope_per_1000_AA / 1000 * separation
        slope_scenarios.append(dict(hypothetical_velocity_slope_m_s_per_1000_AA=slope_per_1000_AA,
            pair_effect_m_s=pair_effect.tolist(), minimum_pair_effect_m_s=float(pair_effect.min()),
            maximum_pair_effect_m_s=float(pair_effect.max()),
            interpretation="Magnitude propagation only; these are hypothetical scenarios, not measured slopes, uncertainty estimates or adopted priors."))
    result = dict(
        interpretation="FIXED-DESIGN DIAGNOSTIC ONLY; no observed offsets or physical amplitude are inferred.",
        background=dict(H0_km_s_Mpc=H0, Omega_m=OM, t0_Gyr=float(age_gyr(0)), flat=True,
                        radiation=False, tstar_s=TSTAR, tstar_fixed_not_fitted=True),
        observable=dict(anchor=2374, comparison=2382, c_m_s=C_M_S,
            reference_wavelengths_AA=wav.tolist(), reference="Oscillator-strength-weighted composite from the same isotope table used by profile fits; no standalone laboratory covariance supplied.",
            D="c ln[(lambda_obs_2382/lambda_obs_2374)/(lambda_ref_2382/lambda_ref_2374)] = v_2382-v_2374",
            Y="ln[(nu_abs_2382/nu_abs_2374)/(nu_lab_2382/nu_lab_2374)] = -D/c, conditional on differential calibration/propagation/gas effects",
            scope="One ratio of two transition-energy differences sharing a lower level, not a complete atomic energy-level distribution.",
            line_center_caveat="Fitted common-velocity offsets are conditional parameters of isotope/Voigt mixtures, not model-free observed centroids."),
        design_targets=rows,
        primary_phenomenological_model="D_i = b + A G_2(z_i) + k_i DeltaLambda_obs_i + e_i + noise_i; G_2(0)=0, G_2(3)=1",
        null="Constant pair offset b; this permits a common laboratory/instrumental offset and is not absolute D=0.",
        physical_identifiability="A is exactly nonidentifiable with unrestricted per-target k_i or e_i. Finite conditional forecasts set those differential terms to known zero; they are not measured physical constraints.",
        practical_calibration_balance="Full-rank unconstrained nuisances prove structural nonidentifiability, not that real calibration errors are arbitrarily large. The close pair suppresses smooth wavelength-slope errors, so finite external calibration constraints can preserve useful conditional precision. No such calibrated per-target prior is measured here.",
        hypothetical_slope_scenarios=slope_scenarios,
        shared_slope_mimic_of_hypothetical_A_100_m_s=dict(
            fitted_common_offset_m_s=float(shared_slope_fit[0]),
            fitted_shared_velocity_slope_m_s_per_1000_AA=float(shared_slope_fit[1] * 1000),
            residual_m_s=shared_slope_residual.tolist(),
            residual_rms_m_s=float(np.sqrt(np.mean(shared_slope_residual**2))),
            residual_max_abs_m_s=float(np.max(np.abs(shared_slope_residual))),
            interpretation="Best unweighted affine approximation to an arbitrary A=100 m/s signal, not a measured instrumental distortion. Its required magnitude must be compared with independent calibration data, not inferred from rank alone."),
        nuisance_projection=rank_results,
        competitor_shapes=comparison,
        equal_weight_forecasts=forecasts,
        common_laboratory_error="A shared laboratory pair error adds sigma_lab^2 11^T; it is absorbed by free b and does not average down. It cancels from cross-target differences and does not degrade slope in this idealized covariance model.",
        sample_span=dict(z_min=float(z.min()), z_max=float(z.max()),
                         age_min_Gyr=float(age_gyr(z.max())), age_max_Gyr=float(age_gyr(z.min())),
                         G2_span=float(np.ptp(g))),
        near_same_age_control=dict(targets=["J053007-250329", "J233156-090802"],
            redshifts=[2.141, 2.143],
            age_separation_Myr=float((age_gyr(2.141) - age_gyr(2.143)) * 1000),
            G2_difference=float(curve(2.143, "inverse_log_p2") - curve(2.141, "inverse_log_p2")),
            expected_D_difference_if_hypothetical_A_100_m_s=float(100 * (curve(2.143, "inverse_log_p2") - curve(2.141, "inverse_log_p2"))),
            interpretation="Under a universal deterministic age-only pair law with the same response, the pair should have almost the same offset. Differences probe the combined sightline, environment, gas-model and calibration terms, not one specific cause; similarities do not establish cosmic variation."),
        no_data_dependent_amplitude_selection=True,
        source_hashes={str(SOURCE.relative_to(ROOT)): sha(SOURCE), str(ATOMIC.relative_to(ROOT)): sha(ATOMIC),
                       str(Path(__file__).relative_to(ROOT)): sha(__file__)},
        software=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))
    (OUT / "design.json").write_text(json.dumps(result, indent=2) + "\n")
    with (OUT / "target_geometry.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    np.savez_compressed(OUT / "design_arrays.npz", redshift=z, G2=g, observed_separation_AA=separation,
                        **{key: value for key, value in designs.items()})
    figure(z, g, result)
    return result


def figure(z, g, result):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    dense = np.linspace(z.min() - .07, z.max() + .07, 250)
    styles = [("inverse_log_p2", "Inverse-log squared", "#202c59", "-"),
              ("inverse_log_p1", "Inverse-log, p=1", "#c65d2d", "--"),
              ("inverse_log_p3", "Inverse-log, p=3", "#20928e", ":"),
              ("log_cosmic_age", "Log cosmic age", "#79579c", "-."),
              ("linear_cosmic_age", "Linear cosmic age", "#767676", "--")]
    for key, label, color, ls in styles:
        axes[0].plot(dense, curve(dense, key), label=label, color=color, ls=ls, lw=1.8)
    axes[0].scatter(z, g, color="#202c59", s=27, zorder=5)
    axes[0].set(xlabel="Absorber redshift", ylabel="Time basis G(z), normalized at z=3", title="(a) The six-target age geometry")
    axes[0].legend(fontsize=8, loc="upper left", frameon=False)
    color = ["#c65d2d", "#20928e", "#79579c", "#767676", "#d99b21"]
    labels = ["p=1", "p=3", "log age", "linear age", "linear z"]
    values = [result["competitor_shapes"][k]["rms_residual_for_hypothetical_A_100_m_s"] for k in result["competitor_shapes"]]
    axes[1].bar(labels, values, color=color)
    axes[1].set_yscale("log")
    axes[1].set(ylabel="RMS unmatched signal (m/s)", title="(b) Curvature left after offset + amplitude")
    axes[1].tick_params(axis="x", rotation=30)
    for j, value in enumerate(values):
        axes[1].annotate(f"{value:.3g}", (j, value), xytext=(0, 4), textcoords="offset points", ha="center", fontsize=8)
    axes[1].text(.02, .97, "Hypothetical A=100 m/s\nEqual weights; no measured signal", transform=axes[1].transAxes,
                 va="top", fontsize=9)
    axes[1].set_ylim(min(values) * .5, max(values) * 25)
    keys = list(result["nuisance_projection"])
    fraction = [100 * result["nuisance_projection"][k]["information_fraction_vs_common_offset"] for k in keys]
    labels = ["Common\noffset", "+ shared\nslope", "Per-target\nslopes", "Per-target\ngas offsets"]
    axes[2].bar(labels, fraction, color=["#202c59", "#20928e", "#b9bec7", "#b9bec7"])
    for j, value in enumerate(fraction):
        axes[2].text(j, value + 2, f"{value:.3g}%" if value else "exactly 0", ha="center", fontsize=9)
    axes[2].set(ylabel="Remaining geometric information (%)", ylim=(0, 120),
                title="(c) Calibration / environment degeneracy")
    axes[2].text(.03, .97, "Unrestricted nuisance amplitudes;\nnot measured calibration errors", transform=axes[2].transAxes, va="top", fontsize=9)
    fig.suptitle("Common Fe II 2374 / 2382 pair: fixed-design diagnostics, not an observed cosmic-time fit", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, .94))
    for ext in ["pdf", "png"]:
        fig.savefig(OUT / f"cross_age_identifiability.{ext}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    summary = generate()
    print(json.dumps(dict(span=summary["sample_span"], nuisance=summary["nuisance_projection"],
                         competitor_shapes=summary["competitor_shapes"], forecasts=summary["equal_weight_forecasts"]), indent=2))
