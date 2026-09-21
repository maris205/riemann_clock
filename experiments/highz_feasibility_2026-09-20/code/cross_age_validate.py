#!/usr/bin/env python3
"""Independent mathematical audit of fixed cross-age design products.

This validator does not import or execute cross_age_design.py. Cosmic ages are
computed by numerical Friedmann-time integration, the competing curve fits use
closed-form centered regression, and forecasts use analytical contrast results.
No observed profile offset is read or interpreted as a physical measurement.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/cross_age"
REPORT = ROOT / "reports/cross_age_independent_review.md"
C = 299792458.0
H0_S = 67.4 * 1000.0 / 3.0856775814913673e22
OM = 0.315
YEAR_S = 365.25 * 86400.0
TP_S = 5.391247e-44
checks = []


def check(name, passed, detail=None):
    row = {"name": name, "pass": bool(passed)}
    if detail is not None:
        row["detail"] = detail
    checks.append(row)


def near(name, observed, expected, atol=2e-11, rtol=2e-10):
    observed, expected = np.asarray(observed), np.asarray(expected)
    check(name, observed.shape == expected.shape and
          np.allclose(observed, expected, atol=atol, rtol=rtol),
          {"maximum_absolute_error": float(np.max(np.abs(observed - expected)))})


def age_seconds(z):
    # dt = da / [a H(a)], integrated from the Big Bang to a=1/(1+z).
    a = 1.0 / (1.0 + z)
    integral, _ = quad(lambda u: np.sqrt(u) / np.sqrt(OM + (1 - OM) * u**3),
                       0.0, a, epsabs=2e-13, epsrel=2e-13)
    return integral / H0_S


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    data = json.loads((OUT / "design.json").read_text())
    rows = data["design_targets"]
    z = np.array([r["z_abs"] for r in rows])
    ages = np.array([age_seconds(v) for v in z])
    t0, t3 = age_seconds(0.0), age_seconds(3.0)
    log0, log3 = np.log(t0 / TP_S), np.log(t3 / TP_S)
    logarithms = np.log(ages / TP_S)
    # Direct ratios independently evaluate the stable production expression.
    def gp(p):
        return ((log0 / logarithms)**p - 1) / ((log0 / log3)**p - 1)
    g = gp(2)
    gc = g - np.mean(g)
    information = float(gc @ gc)
    check("six declared redshifts", np.array_equal(z, [1.629, 1.836, 2.141, 2.143, 2.248, 2.659]))
    check("fixed-design warning", "no observed offsets" in data["interpretation"])
    check("no data-dependent amplitude", data["no_data_dependent_amplitude_selection"] is True)
    near("t0 via Friedmann quadrature", data["background"]["t0_Gyr"], t0 / (1e9 * YEAR_S))
    check("fixed radiation-free background", data["background"]["radiation"] is False and data["background"]["flat"] is True)
    near("declared Planck time", data["background"]["tstar_s"], TP_S, atol=0)
    near("six ages via Friedmann quadrature", [r["cosmic_age_Gyr"] for r in rows], ages / (1e9 * YEAR_S))
    near("six normalized G2 values", [r["G2_z3_normalized"] for r in rows], g)
    near("six unnormalized excesses", [r["h2_minus_1"] for r in rows], (log0 / logarithms)**2 - 1)
    near("G2 endpoint span", data["sample_span"]["G2_span"], np.ptp(g))

    catalogue = json.loads((ROOT / "results/archive_expansion/quality_summary.json").read_text())
    selected = sorted(catalogue["readiness_candidates"], key=lambda r: r["z_abs"])
    check("targets match pre-existing readiness list", [r["target"] for r in rows] == [r["target"] for r in selected])
    near("redshifts match readiness list", z, [r["z_abs"] for r in selected])

    # Parse atomic data directly, without importing any fitting implementation.
    refs = []
    atomic = ROOT / "data/raw/espresso_null/MM_VPFIT_2013-11-10.dat"
    for key in (2374, 2382):
        components = []
        for line in atomic.read_text().splitlines():
            columns = line.split()
            if len(columns) > 5 and columns[0] == "FeII" and key <= float(columns[1]) < key + 1:
                components.append((float(columns[1]), float(columns[2])))
        check(f"four isotope components {key}", len(components) == 4)
        refs.append(sum(w * f for w, f in components) / sum(f for w, f in components))
    refs = np.array(refs)
    delta = (refs[1] - refs[0]) * (1 + z)
    near("laboratory reference wavelengths", data["observable"]["reference_wavelengths_AA"], refs)
    near("observed pair separations", [r["observed_pair_separation_AA"] for r in rows], delta)
    near("observed short wavelengths", [r["observed_2374_AA"] for r in rows], refs[0] * (1 + z))
    near("observed long wavelengths", [r["observed_2382_AA"] for r in rows], refs[1] * (1 + z))
    near("illustrative exact mimic slopes", [r["velocity_slope_m_s_per_AA_that_exactly_mimics_A_100_m_s"] for r in rows], 100 * g / delta)

    # Construct artificial spectral frequencies to verify sign and cancellation.
    # This is an algebraic test, not a mock observation or an injection forecast.
    fractional_log_frequency = np.array([2e-7, -3e-7])
    expected_d = -C * np.diff(fractional_log_frequency)[0]
    shifts, lab_increments = [], []
    lab_errors = np.array([-1e-8, 2e-8])
    for redshift in z:
        observed_lambda = (1 + redshift) * refs * np.exp(-fractional_log_frequency)
        vi = C * np.log(observed_lambda / refs)
        d = vi[1] - vi[0]
        shifts.append(d)
        vwrong = C * np.log(observed_lambda / (refs * np.exp(lab_errors)))
        lab_increments.append((vwrong[1] - vwrong[0]) - d)
    near("positive D for reduced long/short frequency ratio", shifts, np.full(6, expected_d), atol=2e-7)
    near("laboratory pair error is the same for every target", lab_increments, np.full(6, -C * np.diff(lab_errors)[0]), atol=2e-7)
    near("100 m/s fractional frequency-ratio conversion", np.expm1(-100 / C), -3.335640395656554e-7, atol=1e-18)

    # Nuisance rank follows from explicit algebra, without production SVD.
    ranks = data["nuisance_projection"]
    centered_delta = delta - delta.mean()
    shared_slope_residual = gc - centered_delta * (centered_delta @ gc) / (centered_delta @ centered_delta)
    expected_ranks = {
        "common_offset": (1, 2, gc),
        "common_offset_and_shared_wavelength_slope": (2, 3, shared_slope_residual),
        "one_wavelength_slope_per_target": (6, 6, np.zeros(6)),
        "one_environmental_pair_offset_per_target": (6, 6, np.zeros(6))}
    for name, (nr, ar, residual) in expected_ranks.items():
        row = ranks[name]
        check(name + " analytical ranks", row["nuisance_rank"] == nr and row["augmented_rank"] == ar)
        check(name + " additional rank", row["identifiable_additional_rank"] == ar - nr)
        check(name + " structural status", row["structural_nonidentifiability"] == (nr == ar))
        near(name + " response", row["projected_response"], residual)
        near(name + " squared response", row["squared_projected_response"], residual @ residual)
        near(name + " information fraction", row["information_fraction_vs_common_offset"], (residual @ residual) / information)
    for amplitude in (-231.0, 1.0, 100.0, 1000.0):
        slopes = -amplitude * g / delta
        near(f"exact per-target slope cancellation A={amplitude}", amplitude * g + slopes * delta, np.zeros(6), atol=3e-13)

    check("rank result does not assert actual large errors", "not that real calibration errors are arbitrarily large" in data["practical_calibration_balance"])
    for row in data["hypothetical_slope_scenarios"]:
        slope = row["hypothetical_velocity_slope_m_s_per_1000_AA"]
        effects = delta * slope / 1000.0
        near(f"hypothetical slope {slope} propagation", row["pair_effect_m_s"], effects)
        near(f"hypothetical slope {slope} extrema", [row["minimum_pair_effect_m_s"], row["maximum_pair_effect_m_s"]], [effects.min(), effects.max()])
        check(f"hypothetical slope {slope} is not an adopted prior", "not measured slopes, uncertainty estimates or adopted priors" in row["interpretation"])
    affine = data["shared_slope_mimic_of_hypothetical_A_100_m_s"]
    affine_slope = 100 * (centered_delta @ gc) / (centered_delta @ centered_delta)
    affine_intercept = 100 * g.mean() - affine_slope * delta.mean()
    near("shared slope illustrative fitted parameters", [affine["fitted_common_offset_m_s"], affine["fitted_shared_velocity_slope_m_s_per_1000_AA"]], [affine_intercept, affine_slope * 1000], atol=2e-9)
    near("shared slope illustrative residual", affine["residual_m_s"], 100 * shared_slope_residual)
    near("shared slope illustrative RMS", affine["residual_rms_m_s"], 100 * np.sqrt(np.mean(shared_slope_residual**2)))
    near("shared slope illustrative maximum", affine["residual_max_abs_m_s"], 100 * np.max(np.abs(shared_slope_residual)))
    check("shared slope illustration requires external calibration", "not a measured instrumental distortion" in affine["interpretation"])

    alternatives = {
        "inverse_log_p1": gp(1), "inverse_log_p3": gp(3),
        "log_cosmic_age": np.log(t0 / ages) / np.log(t0 / t3),
        "linear_cosmic_age": (t0 - ages) / (t0 - t3),
        "linear_redshift": z / 3.0}
    for name, alternative in alternatives.items():
        row = data["competitor_shapes"][name]
        ac = alternative - alternative.mean()
        slope = (ac @ gc) / (ac @ ac)
        intercept = g.mean() - slope * alternative.mean()
        residual = gc - slope * ac
        norm = float(np.linalg.norm(residual))
        near(name + " regression parameters", [row["fitted_offset"], row["fitted_amplitude"]], [intercept, slope])
        near(name + " residual", row["residual_after_offset_and_amplitude"], residual)
        near(name + " norm", row["residual_norm"], norm)
        near(name + " centered correlation", row["centered_correlation"], (gc @ ac) / np.linalg.norm(gc) / np.linalg.norm(ac))
        near(name + " relative curvature", row["fractional_shape_residual_relative_to_centered_p2"], norm / np.linalg.norm(gc))
        near(name + " arbitrary100m/s maximum", row["max_abs_residual_for_hypothetical_A_100_m_s"], 100 * np.max(np.abs(residual)))
        near(name + " arbitrary100m/s RMS", row["rms_residual_for_hypothetical_A_100_m_s"], 100 * norm / np.sqrt(6))
        sig = row["hypothetical_sigma_each_for_expected_delta_chi2_9_at_A_100_m_s"]
        near(name + " Gaussian mean separation only", (100 * residual / sig) @ (100 * residual / sig), 9, atol=2e-7)
        check(name + " no observed-amplitude/significance claim", "arbitrary, not fitted or predicted" in row["warning"] and "No detection significance" in row["warning"])

    # With free intercept, a common rank-one covariance affects the intercept,
    # whereas the centered slope contrast has variance sigma^2 / sum(g-gbar)^2.
    for j, row in enumerate(data["equal_weight_forecasts"]):
        sigma_a = row["sigma_independent_pair_m_s"] / np.sqrt(information)
        near(f"forecast {j} slope analytical contrast", row["conditional_sigma_A_m_s"], sigma_a)
        near(f"forecast {j} endpoint analytical contrast", row["conditional_sigma_sample_endpoint_difference_m_s"], sigma_a * np.ptp(g))
        common = row["sigma_shared_laboratory_pair_m_s"]**2 * np.ones((6, 6))
        near(f"forecast {j} common lab cancels contrast", gc @ common @ gc, 0.0, atol=1e-12)

    pair = data["near_same_age_control"]
    ix = [next(i for i, r in enumerate(rows) if r["target"] == target) for target in pair["targets"]]
    near("same-age control redshifts", pair["redshifts"], z[ix])
    near("same-age control age difference", pair["age_separation_Myr"], abs(np.diff(ages[ix])[0]) / (YEAR_S * 1e6))
    near("same-age control response difference", pair["G2_difference"], abs(np.diff(g[ix])[0]))
    near("same-age control illustrative velocity difference", pair["expected_D_difference_if_hypothetical_A_100_m_s"], 100 * abs(np.diff(g[ix])[0]))
    check("same-age control does not establish variation", "similarities do not establish cosmic variation" in pair["interpretation"])

    with (OUT / "target_geometry.csv").open(newline="") as handle:
        csv_rows = list(csv.DictReader(handle))
    check("CSV target order", [r["target"] for r in csv_rows] == [r["target"] for r in rows])
    for key in rows[0]:
        if key != "target":
            near("CSV " + key, [float(r[key]) for r in csv_rows], [r[key] for r in rows])
    with np.load(OUT / "design_arrays.npz") as arrays:
        near("NPZ redshifts", arrays["redshift"], z)
        near("NPZ time basis", arrays["G2"], g)
        near("NPZ wavelength separation", arrays["observed_separation_AA"], delta)
        near("NPZ common offset", arrays["common_offset"], np.ones((6, 1)))
        near("NPZ common slope", arrays["common_offset_and_shared_wavelength_slope"], np.column_stack([np.ones(6), delta]))
        near("NPZ per-target slopes", arrays["one_wavelength_slope_per_target"], np.column_stack([np.ones(6), np.diag(delta)]))
        near("NPZ per-target gas offsets", arrays["one_environmental_pair_offset_per_target"], np.column_stack([np.ones(6), np.eye(6)]))
    for filename, expected in data["source_hashes"].items():
        check("source hash " + filename, digest(ROOT / filename) == expected)
    for extension in ("png", "pdf"):
        check("figure exists " + extension, (OUT / f"cross_age_identifiability.{extension}").stat().st_size > 10000)

    failed = [row for row in checks if not row["pass"]]
    result = dict(status="PASS" if not failed else "FAIL", checks_passed=len(checks) - len(failed),
                  checks_total=len(checks), failed=failed, checks=checks,
                  scope="Independent fixed-design mathematics and artifact audit only. No measured cosmic-time amplitude or discovery significance is inferred.",
                  independent_methods=["Cosmic age by numerical Friedmann-time integration", "Closed-form centered regressions", "Explicit nuisance reparameterization", "Analytical Gaussian contrast variance", "Constructed wavelength/frequency sign and shared-laboratory-error check"],
                  source_hashes={"code/cross_age_design.py": digest(ROOT / "code/cross_age_design.py"),
                                 "code/cross_age_validate.py": digest(Path(__file__)),
                                 "results/cross_age/design.json": digest(OUT / "design.json")},
                  python=platform.python_version())
    (OUT / "independent_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    REPORT.write_text(f"""# Cross-age fixed-design independent review

**{result['status']}: {result['checks_passed']}/{result['checks_total']} numerical and artifact checks.**

This review concerns the geometry of six pre-existing candidate redshifts. It does not fit observed line offsets, establish a physical amplitude, or assign discovery significance. The numerical reviewer did not import or execute the production generator. Ages were independently obtained by numerical Friedmann-time integration; curve comparisons used closed-form centered regression; nuisance degeneracies and Gaussian contrast errors were checked analytically.

The sign convention is consistent: for the 2382-minus-2374 log-wavelength velocity contrast, $D=-c\\ln[(\\nu_{{2382}}/\\nu_{{2374}})_z/(\\nu_{{2382}}/\\nu_{{2374}})_{{\\rm lab}}]$. Positive $D$ means a smaller frequency ratio. Photon transition-energy differences obey the same ratio; an entire atomic level distribution is not measured. At an illustrative +100 m/s, the relative frequency-ratio change is -3.3356404e-7. Laboratory wavelength errors shared by the targets give one common pair offset, not independent target errors. Such a common term cancels from target contrasts and is absorbed by a free intercept in the idealized covariance forecast.

Independent ages span {ages.min() / (1e9 * YEAR_S):.8f} to {ages.max() / (1e9 * YEAR_S):.8f} Gyr for the declared fixed flat matter-plus-Lambda background. Radiation is neglected. The Planck-time reference is a design assumption, not a derived atomic response.

The per-target slope design is exactly full row rank: $D_i=b+A G_2(z_i)+s_i\\Delta\\lambda_i$ is unchanged under $A\\mapsto A+\\delta A$, $s_i\\mapsto s_i-\\delta A G_2(z_i)/\\Delta\\lambda_i$. Every pair has nonzero wavelength separation. Consequently, an unrestricted slope or environmental pair offset for each target removes physical identifiability. Assuming one common slope is much more restrictive and retains only {100 * ranks['common_offset_and_shared_wavelength_slope']['information_fraction_vs_common_offset']:.6f}% of the geometric response information relative to a common offset alone. Finite forecast errors condition on differential nuisance terms being known, and must not be presented as measured constraints.

This rank statement does not establish that actual calibration errors have arbitrary magnitude. The close pair mitigates a smooth wavelength slope: illustrative slopes of 50, 200 and 500 m/s per 1000 Angstrom give pair biases of 1.09--1.52, 4.37--6.08 and 10.92--15.19 m/s, respectively. These are propagation scenarios, not measured errors or adopted priors. Fitting an intercept and common slope to the arbitrary A=100 m/s curve requires {affine['fitted_shared_velocity_slope_m_s_per_1000_AA']:.6f} m/s per 1000 Angstrom and leaves {affine['residual_rms_m_s']:.8f} m/s RMS mismatch. A=100 m/s denotes the response at z=3 relative to z=0; its difference across these six targets is only {100 * np.ptp(g):.8f} m/s. Actual external calibration information is required to assess whether any relevant distortion is plausible.

After refitting both offset and amplitude, the p=1 and p=3 inverse-log competitors leave only {data['competitor_shapes']['inverse_log_p1']['rms_residual_for_hypothetical_A_100_m_s']:.8f} and {data['competitor_shapes']['inverse_log_p3']['rms_residual_for_hypothetical_A_100_m_s']:.8f} m/s RMS mismatch for the arbitrary normalization A=100 m/s. These are noiseless, equal-weight design comparisons, not physical predictions, likelihood evidence, a statistical power calculation, or a significance estimate.

The near-equal-age sightlines J053007-250329 (z=2.141) and J233156-090802 (z=2.143) differ by {pair['age_separation_Myr']:.8f} Myr. A universal deterministic age-only response with A=100 m/s would differ by only {pair['expected_D_difference_if_hypothetical_A_100_m_s']:.8f} m/s. Comparing these targets is therefore an informative combined environment, calibration and gas-model control. Agreement would not establish cosmic variation; disagreement would not uniquely identify which nuisance term is responsible.

Reference wavelengths, source hashes, CSV/NPZ arrays, age conversion, response curves, covariance calculations and the existence of figure outputs were checked. One initially unreachable optional curve branch was reported and fixed before the audited products were finalized; it had not entered the generated comparisons. The two transitions' common lower state is independently supported by [Bainbridge and Webb (2017), Table A1](https://academic.oup.com/view-large/64632789). The full check list is in `results/cross_age/independent_validation.json`.
""")
    print(json.dumps({key: result[key] for key in ("status", "checks_passed", "checks_total", "failed")}, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
