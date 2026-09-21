#!/usr/bin/env python3
"""Independent checks of the conditional forecast; never imports the forecast code.

Independent implementations:
- cosmic age from quadrature in scale factor rather than the asinh formula;
- 60-digit direct logarithm powers rather than expm1/log1p;
- ordinary least-squares projection rather than QR;
- regression of six independent stratum means rather than N x N GLS;
- line-level Monte Carlo with explicit nuisance shifts/slopes and group errors.

The forecast is a toy design with unvalidated atomic response coefficients.
Passing these checks does not validate a physical Riemann-to-spectrum mapping.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import platform

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm, beta

BASE = Path(__file__).resolve().parents[1]
SUMMARY = BASE / "results/forecast_summary.json"
OUT_JSON = BASE / "results/forecast_validation.json"
OUT_MD = BASE / "reports/forecast_validation.md"
MPC_KM = 3.0856775814913673e19
YEAR = 31557600.0
SEED = 71620493
NMC = 10000
mp.mp.dps = 60
checks = []


def record(name, success, detail):
    checks.append({"check": name, "pass": bool(success), "detail": detail})


def age_quadrature(z, h0, om):
    # t(a)=H0^-1 integral_0^a sqrt(a')/sqrt(Om+(1-Om)*a'^3) da'.
    a = 1.0 / (1.0 + float(z))
    integ, err = quad(lambda aa: np.sqrt(aa) / np.sqrt(om+(1-om)*aa**3),
                      0.0, a, epsabs=1e-13, epsrel=3e-13)
    return integ * MPC_KM / h0 / YEAR / 1e9


def hp_direct(z, p, h0, om, tstar):
    t0 = mp.mpf(str(age_quadrature(0, h0, om))) * mp.mpf(str(YEAR)) * 10**9
    tz = mp.mpf(str(age_quadrature(z, h0, om))) * mp.mpf(str(YEAR)) * 10**9
    return (mp.log(t0/mp.mpf(str(tstar)))/mp.log(tz/mp.mpf(str(tstar))))**p - 1


def normalized(z, p, h0, om, tstar):
    denominator = hp_direct(3, p, h0, om, tstar)
    return np.array([float(hp_direct(zi, p, h0, om, tstar)/denominator) for zi in np.atleast_1d(z)])


def orthogonal(k, wave, degree):
    # Different centering/scaling is harmless: the polynomial span is unchanged.
    xx = np.log(wave / wave[0])
    xx /= np.max(xx)
    design = np.vander(xx, degree+1, increasing=True)
    fitted = design @ np.linalg.lstsq(design, k, rcond=None)[0]
    return k-fitted, xx, design


def interval_binomial(k, n, confidence=.999):
    alpha = 1-confidence
    return [0. if k == 0 else float(beta.ppf(alpha/2, k, n-k+1)),
            1. if k == n else float(beta.ppf(1-alpha/2, k+1, n-k))]


def run():
    summary = json.loads(SUMMARY.read_text())
    h0 = summary['background']['H0_km_s_Mpc']
    om = summary['background']['Omega_m']
    ts = summary['background']['tstar_s']
    k = np.array(summary['synthetic_response']['K'])
    wave = np.array(summary['synthetic_response']['rest_vacuum_wavelength_grid_A_lookup_only'])
    record('RMS response normalization', np.isclose(np.mean(k*k), 1, rtol=2e-14), float(np.mean(k*k)))
    ages = []
    for item in summary['time_leverage']:
        z = item['z']
        independent_age = age_quadrature(z, h0, om)
        direct_excess = float(hp_direct(z, 2, h0, om, ts))
        direct_shape = float(normalized([z], 2, h0, om, ts)[0])
        agree = (np.isclose(independent_age, item['age_Gyr'], rtol=2e-12) and
                 np.isclose(direct_excess, item['h2_minus_one'], rtol=2e-12, atol=1e-14) and
                 np.isclose(direct_shape, item['G2_z3_normalized'], rtol=2e-12, atol=1e-14))
        ages.append({'z': z, 'quadrature_age_Gyr': independent_age,
                     'direct_high_precision_h2_minus_one': direct_excess,
                     'G2_normalized': direct_shape,
                     'age_difference_Gyr': independent_age-item['age_Gyr']})
        record(f'age and direct log powers z={z}', agree, ages[-1])
    projection = {}
    for deg in (1, 2):
        perp, xx, nuisance = orthogonal(k, wave, deg)
        s2 = float(perp@perp)
        expected = summary['response_retained_squared_norm'][str(deg)]
        projection[str(deg)] = {'norm2': s2, 'fraction': s2/float(k@k),
                                'max_nuisance_inner_product': float(np.max(np.abs(nuisance.T@perp)))}
        record(f'independent polynomial projection degree {deg}',
               np.isclose(s2, expected['norm2'], rtol=2e-13) and
               np.allclose(perp, expected['projected_K'], atol=3e-13, rtol=2e-12), projection[str(deg)])
    for label, response in [('common_shift', np.ones(6)), ('log_wavelength_slope', np.log(wave))]:
        pp, _, _ = orthogonal(response, wave, 1)
        record(f'exact degeneracy {label}', np.linalg.norm(pp) < 1e-12, {'residual_norm': float(np.linalg.norm(pp))})
    perp, xx, nuisance = orthogonal(k, wave, 1)
    norm2 = float(perp@perp)
    rng = np.random.default_rng(SEED)
    output_scenarios = {}
    for name, item in summary['scenarios'].items():
        n = item['n_independent_sightlines_assumed']
        m = n//6
        z6 = np.linspace(*item['z_range'], 6)
        g6 = normalized(z6, 2, h0, om, ts)
        ss = item['sigma_line_mps_assumed']
        floor = item['six_bin_correlated_floor_mps_assumed']
        variance_bin = ss**2/norm2/m + floor**2
        centered = g6-np.mean(g6)
        slope_weights = centered/(centered@centered)
        sigma_b = np.sqrt(variance_bin/(centered@centered))
        expected_b = item['sigma_B_at_z3_mps']
        ci95 = norm.ppf(.975)*sigma_b
        target_snr = brentq(lambda snr: norm.sf(5-snr)+norm.cdf(-5-snr)-.9, 5, 10, xtol=1e-13)
        power90_b = target_snr*sigma_b
        scalar_checks = {
            'sigma_B': bool(np.isclose(sigma_b, expected_b, rtol=3e-12)),
            'null_95_halfwidth': bool(np.isclose(ci95, item['null_95_interval_halfwidth_B_mps'], rtol=3e-12)),
            '5sigma90power': bool(np.isclose(power90_b, item['B_mps_for_5sigma_test_90percent_power'], rtol=3e-12)),
        }
        record(f'six-bin sufficient-statistic regression {name}', all(scalar_checks.values()), scalar_checks)
        # An arbitrary unknown offset for every stratum absorbs the signal completely.
        memberships = np.repeat(np.eye(6), m, axis=0)
        gn = np.repeat(g6, m)
        unrestricted_rank = int(np.linalg.matrix_rank(memberships))
        appended_rank = int(np.linalg.matrix_rank(np.column_stack([memberships, gn])))
        record(f'unconstrained bin offsets destroy identifiability {name}', unrestricted_rank == appended_rank == 6,
               {'bin_nuisance_rank': unrestricted_rank, 'with_cosmic_shape_rank': appended_rank})
        # Check trend-aligned covariance treatment independently at group-mean level.
        group_cov = np.eye(6)*variance_bin + floor**2*np.outer(g6, g6)
        x6 = np.column_stack([np.ones(6), g6])
        params_cov = np.linalg.inv(x6.T@np.linalg.solve(group_cov, x6))
        sigma_align = np.sqrt(params_cov[1, 1])
        record(f'extra aligned random floor {name}',
               np.isclose(sigma_align, item['sigma_B_with_extra_trend_aligned_floor_mps'], rtol=3e-12),
               {'sigma_B_aligned_direct_covariance': float(sigma_align),
                'interpretation': 'Known zero-mean Gaussian random-effect variance, not a bound for an arbitrary unknown systematic.'})
        shape_checks = {}
        for p in (1, 3):
            gp = normalized(z6, p, h0, om, ts)
            design_p = np.column_stack([np.ones(6), gp])
            resid = g6-design_p@np.linalg.lstsq(design_p, g6, rcond=None)[0]
            threshold = 5*np.sqrt(variance_bin/(resid@resid))
            expected = item['signal_shape_discrimination'][str(p)]
            rms = np.sqrt(np.mean(resid**2))
            shape_checks[str(p)] = {'B_for_expected_delta_chi2_25_mps': float(threshold),
                                    'rms_after_offset_and_amplitude': float(rms)}
            record(f'shape p={p} conditional distance {name}',
                   np.isclose(threshold, expected['B_mps_for_expected_delta_chi2_25'], rtol=3e-8) and
                   np.isclose(rms, expected['rms_shape_residual_after_offset_and_amplitude'], rtol=3e-8), shape_checks[str(p)])
        # Reconstruct independent full line data rather than sampling fitted B.
        # Each line has an independent error; all objects in a stratum share one mode error.
        # Large arbitrary line-common and log-wave nuisance offsets must cancel.
        simulations = []
        for snr in (0., target_snr):
            truth = snr*sigma_b
            errors = rng.normal(0, ss, size=(NMC, n, 6))
            block_errors = np.repeat(rng.normal(0, floor, size=(NMC, 6)), m, axis=1)
            redshift_nuisance = rng.normal(0, 1e4, size=(NMC, n, 1))
            slope_nuisance = rng.normal(0, 1e3, size=(NMC, n, 1))
            lines = ((truth*gn[None, :, None]+block_errors[:, :, None])*k[None, None, :]
                     + redshift_nuisance + slope_nuisance*xx[None, None, :] + errors)
            mode_estimates = (lines@perp)/norm2
            group_means = mode_estimates.reshape(NMC, 6, m).mean(axis=2)
            b_est = group_means@slope_weights
            residual = (b_est-truth)/sigma_b
            coverage_count = int(np.sum(np.abs(residual) < norm.ppf(.975)))
            rejection_count = int(np.sum(np.abs(b_est/sigma_b)>5))
            coverage_interval = interval_binomial(coverage_count, NMC)
            expected_power = float(norm.sf(5-snr)+norm.cdf(-5-snr))
            power_interval = interval_binomial(rejection_count, NMC)
            this = {'injected_B_mps': float(truth), 'n_simulations': NMC,
                    'residual_mean_sigma': float(residual.mean()),
                    'residual_sd_sigma': float(residual.std(ddof=1)),
                    'coverage_95': coverage_count/NMC, 'coverage_99_9percent_binomial_interval': coverage_interval,
                    'rejection_abs_snr_gt5': rejection_count/NMC,
                    'analytic_rejection_probability': expected_power,
                    'rejection_99_9percent_binomial_interval': power_interval}
            simulations.append(this)
            sim_ok = (abs(residual.mean()) < 4/np.sqrt(NMC) and
                      abs(residual.std(ddof=1)-1) < 4/np.sqrt(2*(NMC-1)) and
                      coverage_interval[0] <= .95 <= coverage_interval[1] and
                      power_interval[0] <= expected_power <= power_interval[1])
            record(f'independent line-level MC snr={snr:.6g} {name}', sim_ok, this)
        # Verify original MC results against broad statistical fluctuation bounds.
        for j, original in enumerate(summary['monte_carlo'][name]):
            nn = original['n_simulations']
            coh = (abs(original['mean_bias_in_sigma']) < 4/np.sqrt(nn) and
                   abs(original['residual_std_in_sigma']-1) < 4/np.sqrt(2*(nn-1)) and
                   abs(original['coverage_95']-.95) < 4*np.sqrt(.95*.05/nn))
            if original['injected_snr'] != 0:
                p = original['analytic_fraction_abs_snr_gt5']
                coh = coh and abs(original['fraction_abs_snr_gt5']-p) < 4*np.sqrt(p*(1-p)/nn)
            record(f'original MC consistency entry {j} {name}', coh,
                   {'n_simulations': nn, 'interpretation': 'Sampling consistency, not empirical verification of a 5-sigma null tail.'})
        output_scenarios[name] = {'sigma_B_mps_independent': float(sigma_b),
                    'null_95_halfwidth_B_mps': float(ci95),
                    'five_sigma_90percent_power_B_mps': float(power90_b),
                    'variance_single_bin_mean_mps2': float(variance_bin),
                    'shape_separation': shape_checks, 'line_level_monte_carlo': simulations}
    # Check whether the original uncertainty-floor model obeys saturation rather than 1/sqrt(N).
    floor_limits = {}
    for name, item in summary['scenarios'].items():
        g = normalized(np.linspace(*item['z_range'], 6), 2, h0, om, ts)
        denom = np.sqrt(np.sum((g-g.mean())**2))
        floor_limits[name] = item['six_bin_correlated_floor_mps_assumed']/denom
    limitations = [
        'The response K is invented for an experimental-design calculation. It is not an atomic calculation, an alpha sensitivity, or an observed signal.',
        'The forecast transfers a logarithmic time template to a deterministic line-centroid effect. The paper\'s extra error standard deviation is a different observable, and no equivalence has been derived.',
        'The six redshift-bin floors have known Gaussian variances and zero-mean independent draws. Arbitrary free bin offsets or an arbitrary trend-aligned bias are unidentifiable from the proposed signal.',
        'Equal independent centroid errors, six usable lines per object, equal bin populations, fixed K, and a correctly specified nuisance span are optimistic design assumptions requiring measurement in the pilot.',
        'The polynomial nuisance is in log wavelength. Real intra-order distortions, blends, isotope mixtures, velocity subcomponents, and model-dependent covariance are not removed by this algebra automatically.',
        'The Planck-time tstar, fixed exponent, and matter+Lambda age coordinate are assumptions. Normalizing at z=3 removes amplitude dependence, not physical uncertainty about the response law.',
        'The 5000- or 10000-run Monte Carlo can assess ordinary coverage and approximately 90% power. It cannot establish a 5.7e-7 null false-positive rate; that tail is analytic under the ideal Gaussian model.',
        'Shape separation thresholds use the expected squared residual against one fixed competitor after refitting its amplitude and intercept. They are not model-selection discovery significances.',
        'Shape-separation amplitudes of 69 to 1632 km/s extrapolate the toy centroid model far beyond a validated response domain. They are evidence of weak exponent discrimination, not proposed physical target amplitudes.',
        'The future-instrument scenario is conditional on delivered coverage, sensitivity, suitable objects, and systematic control; its name or numbers do not guarantee a commissioning date or a 2036 result.'
    ]
    report = {'status': 'PASS' if all(c['pass'] for c in checks) else 'FAIL',
              'created_utc': datetime.now(timezone.utc).isoformat(),
              'method': 'Independent quadrature, high-precision direct powers, six-bin regression and line-level injection; no import of forecast_sensitivity.py',
              'source_summary_sha256': hashlib.sha256(SUMMARY.read_bytes()).hexdigest(),
              'source_code_sha256': hashlib.sha256((BASE/'code/forecast_sensitivity.py').read_bytes()).hexdigest(),
              'software': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__, 'mpmath': mp.__version__},
              'monte_carlo_seed': SEED, 'n_checks': len(checks), 'checks': checks,
              'age_and_log_law': ages, 'projections': projection, 'scenarios': output_scenarios,
              'sigma_B_infinite_N_at_fixed_six_bin_floor_mps': floor_limits,
              'interpretive_limitations': limitations}
    OUT_JSON.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    lines = ['# Independent validation of the conditional high-redshift forecast', '',
             f'**Result: {report["status"]}; {len(checks)} numerical checks.** This is validation of the planning arithmetic, not a validation of the physical hypothesis.', '',
             'The validator does not import the forecast program. It computes cosmic age by numerical integration in scale factor; evaluates the logarithmic law with 60-digit direct powers; projects nuisance modes by least squares; collapses the covariance into six independent bin means; and simulates all six line measurements with explicit offsets, calibration slopes and shared bin errors.', '',
             f'At z = 3, the age is {ages[5]["quadrature_age_Gyr"]:.9f} Gyr and h₂ − 1 = {ages[5]["direct_high_precision_h2_minus_one"]:.10f}, or {100*ages[5]["direct_high_precision_h2_minus_one"]:.6f}% of the assumed present-day amplitude. G₂(3) = 1 by definition. This percentage is not a fractional measured change in any physical constant.', '',
             '| Conditional scenario | Independent σ(B), m/s | Null 95% half-width, m/s | B for 5σ threshold at 90% power, m/s | Fixed-six-bin floor as N→∞, m/s |',
             '|---|---:|---:|---:|---:|']
    for name, sc in output_scenarios.items():
        lines.append(f'| {name} | {sc["sigma_B_mps_independent"]:.6f} | {sc["null_95_halfwidth_B_mps"]:.6f} | {sc["five_sigma_90percent_power_B_mps"]:.6f} | {floor_limits[name]:.6f} |')
    lines.extend(['', 'The known-covariance GLS is exactly reproduced by six equal-variance bin means: Var(B) = [σ_line²/(m ‖K⊥‖²) + σ_floor²] / Σ_bin(G − mean(G))². The fitted intercept is essential; it removes any constant response offset. The stated 90% power amplitude uses the two-sided Gaussian rejection probability, solved independently by root finding.', '',
                  f'The response retains {100*projection["1"]["fraction"]:.6f}% of its squared norm after constant plus log-wavelength slope removal, and {100*projection["2"]["fraction"]:.6f}% after adding curvature. A common shift or pure calibration slope has zero information. If all six bin offsets are free, the cosmic-time column lies in their span and B cannot be identified.', '',
                  f'The independent line-level simulation uses {NMC:,} trials per scenario at zero signal and at the 90%-power amplitude. Large random common shifts and wavelength slopes cancel as expected. All reported intervals and dispersions are statistically compatible with their Gaussian design targets. The complete Monte Carlo results and 99.9% binomial intervals are in `results/forecast_validation.json`.', '',
                  '## Scientific qualifications that must accompany use', ''])
    for item in limitations:
        lines.append('- '+item)
    lines.extend(['', '## Reproduction', '', '```bash', 'python code/forecast_sensitivity.py', 'python code/validate_forecast.py', '```', '',
                  f'Input summary SHA-256: `{report["source_summary_sha256"]}`.',
                  f'Forecast code SHA-256: `{report["source_code_sha256"]}`.', '',
                  'No observation was fitted in this calculation. A pass supports the internal numerical calculation conditional on its design assumptions only.', ''])
    if report['status'] != 'PASS':
        lines.extend(['## Failed checks', '']+[f'- {c["check"]}: {c["detail"]}' for c in checks if not c['pass']])
    OUT_MD.write_text('\n'.join(lines))
    print(json.dumps({'status': report['status'], 'checks': len(checks),
                      'failed_checks': [c['check'] for c in checks if not c['pass']],
                      'independent_sigmas_B': {k:v['sigma_B_mps_independent'] for k,v in output_scenarios.items()}}, indent=2))
    if report['status'] != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    run()
