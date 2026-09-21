"""Audit the focused spectroscopy delivery; no nonlinear fits or old-paper rebuild.

This is an integrity/export audit, not a new physical validation, significance
calibration, independent peer review, or assessment of publication acceptance.
Run after exporting tables, building the PDF, and producing its PDF preflight.
"""
from __future__ import annotations
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import math
import re
import sys
from urllib.parse import unquote

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / 'experiments/highz_feasibility_2026-09-20'
PAPER = ROOT / 'paper_spectroscopy'
BUILD = ROOT / 'build_spectroscopy'
REPORT = ROOT / 'reports'
CHECKS = []
INPUTS = {}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(name, passed, detail=None):
    row = {'name': name, 'passed': bool(passed)}
    if detail is not None:
        row['detail'] = detail
    CHECKS.append(row)


def read(path):
    path = Path(path)
    INPUTS[str(path.relative_to(ROOT))] = sha(path)
    return json.loads(path.read_text())


def text(path):
    path = Path(path)
    INPUTS[str(path.relative_to(ROOT))] = sha(path)
    return path.read_text()


def close(x, y, atol=1e-8):
    return math.isclose(float(x), float(y), abs_tol=atol, rel_tol=1e-9)


def allclose(x, y, atol=1e-8):
    return np.allclose(x, y, atol=atol, rtol=1e-9)


def audit(name, fn):
    try:
        fn()
    except Exception as exc:
        record(name + ': completed without missing or malformed dependency', False,
               f'{type(exc).__name__}: {exc}')


def hash_mapping(mapping, base):
    failures = []
    for rel, expected in mapping.items():
        path = base / rel
        if not path.is_file() or sha(path) != expected:
            failures.append(rel)
    return failures


def verify_tex():
    files = sorted(PAPER.glob('*.tex'))
    texts = {p.name: text(p) for p in files}
    combined = '\n'.join(texts.values())
    missing_inputs, missing_figures = [], []
    for content in texts.values():
        for name in re.findall(r'\\input\{([^}]+)\}', content):
            p = PAPER / name
            if not p.suffix:
                p = p.with_suffix('.tex')
            if not p.is_file():
                missing_inputs.append(name)
        for name in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', content):
            candidates = [PAPER / name, PAPER / 'figures' / name, ROOT / 'figures' / name]
            if not any(p.is_file() for p in candidates):
                missing_figures.append(name)
    record('TeX input files resolve', not missing_inputs, missing_inputs)
    record('TeX figure files resolve', not missing_figures, missing_figures)
    bib = text(PAPER / 'references.bib')
    bibkeys = re.findall(r'@\w+\s*\{\s*([^,\s]+)', bib)
    citations = set()
    for group in re.findall(r'\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]+)\}', combined):
        citations.update(k.strip() for k in group.split(','))
    record('Citations resolve and bibliography keys are unique',
           not (citations - set(bibkeys)) and len(bibkeys) == len(set(bibkeys)),
           {'cited_keys': sorted(citations), 'missing': sorted(citations - set(bibkeys))})
    labels = re.findall(r'\\label\{([^}]+)\}', combined)
    refs = re.findall(r'\\(?:eqref|ref|autoref)\{([^}]+)\}', combined)
    record('Source references resolve and labels are unique',
           not (set(refs) - set(labels)) and len(labels) == len(set(labels)),
           {'missing': sorted(set(refs) - set(labels))})
    record('Focused manuscript is a single-column article',
           r'\documentclass[11pt,a4paper]{article}' in texts['main.tex']
           and 'twocolumn' not in combined)
    record('Same-transition order sign is explicit',
           'lower-index minus higher-index' in texts['order_results.tex'])


def verify_controls():
    provenance = read(PAPER / 'table_provenance.json')
    rows = provenance['control_rows']
    expected = []
    b = read(EXP / 'results/espresso_null/comparison.json')['comparisons']
    for case in ['primary', 'expected_fluctuation']:
        x = next(r for r in b if r['case'] == case)
        expected.append((x['n_pixels'], x['null_parameters'], x['null_chi2'], x['delta_chi2']))
    for folder, prefix in [('empirical_gls', ''), ('joint_controls', 'core01_gls_'),
                           ('joint_controls', 'core03_gls_'), ('gas_structure', '')]:
        x = read(EXP / f'results/{folder}/{prefix}comparison.json')
        n = read(EXP / f'results/{folder}/{x["null"]}.json')
        expected.append((n['ndata'], n['npar'], n['chi2'], x['delta_chi2']))
    record('Six control comparisons reproduce saved source numbers',
           len(rows) == 6 and all(allclose(r[1:], e) for r, e in zip(rows, expected)))
    table = text(PAPER / 'control_table.tex')
    absent = []
    for label, n, k, q, dq in rows:
        fragment = f'{label} & {n} & {k} & {q:.2f} & {dq:.2f}'
        if fragment not in table:
            absent.append(fragment)
    record('Control TeX table reproduces exported rounded rows', not absent, absent)
    bad = hash_mapping(provenance['experiment_source_sha256'], EXP)
    record('Table-export input checksums match current sources', not bad, bad)
    record('Table-export code checksum matches',
           provenance['export_code_sha256'] == sha(ROOT / 'code/export_spectroscopy_tables.py'))


def verify_common_pair():
    summary = read(EXP / 'results/common_pair/refined_summary.json')
    harmonized = read(EXP / 'results/common_pair/harmonized_measurements.json')
    hs = {r['target']: r for r in harmonized['comparisons']}
    for row in summary['comparisons']:
        target = row['target']
        h = hs[target]
        cmp = row['full_comparison']
        n = read(EXP / cmp['null_file'])
        a = read(EXP / cmp['alternative_file'])
        record(f'{target}: common-pair full H0/H1 source identity and arithmetic',
               close(n['chi2'], cmp['chi2_null']) and close(a['chi2'], cmp['chi2_alternative'])
               and close(n['chi2'] - a['chi2'], cmp['delta_chi2'])
               and cmp['added_parameters'] == len(row['lines']) - 1
               and close(h['D_m_s'], row['D_m_s'])
               and close(h['conditional_sigma_m_s'], row['conditional_sigma_m_s'])
               and row['unrestricted_reference_gap_chi2'] <= row['baseline_tolerance_chi2'])
        bad = hash_mapping(row['source_sha256'], EXP)
        ref = EXP / row['point_estimate_reference_file']
        record(f'{target}: source and selected-point hashes',
               not bad and sha(ref) == row['point_estimate_reference_sha256'], bad)
        x = np.asarray(row['profile_grid_D_m_s'])
        y = np.asarray(row['profile_delta_chi2'])
        zero = int(np.argmin(y))
        crossing = []
        for side in [range(zero - 1, -1, -1), range(zero, len(x) - 1)]:
            answer = None
            for i in side:
                if (y[i] - 1) * (y[i + 1] - 1) <= 0 and y[i] != y[i + 1]:
                    answer = float(x[i] + (1-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]))
                    break
            crossing.append(answer)
        record(f'{target}: connected delta-Q=1 crossing interpolation',
               all(v is not None for v in crossing)
               and allclose(crossing, h['profile_delta1_interval_m_s']), crossing)
        record(f'{target}: conditional-only interpretation remains encoded',
               not h['precision_measurement'] and not h['physical_time_inference']
               and h['other_relative_line_shifts_float']
               and h['profile_is_conditional_connected_contour'])
    p = read(EXP / 'results/common_pair/export_provenance.json')
    identities = {
        'input_summary_sha256': EXP / 'results/common_pair/refined_summary.json',
        'export_script_sha256': EXP / 'code/common_pair_report.py',
        'harmonized_sha256': EXP / 'results/common_pair/harmonized_measurements.json',
        'report_sha256': EXP / 'reports/common_pair_results_cn.md',
        'figure_pdf_sha256': EXP / 'results/common_pair/common_pair_profiles.pdf'}
    record('Common-pair final export identities match', all(p[k] == sha(v) for k, v in identities.items()))


def verify_order():
    s = read(EXP / 'results/order_response/model_summary.json')
    source = text(PAPER / 'order_results.tex')
    rows_ok, failure = True, []
    for line in ['2600', '2382']:
        for model in ['G0', 'G1', 'A0', 'A1']:
            r = s[line]['full'][model]
            c0 = s[line]['full'][model[0] + '0']
            q = r['chi2']
            expected = f'{line} & {model[0]}$_{model[1]}$ & {q:.3f}'
            ok = expected in source and close(sum(t['chi2'] for t in r['per_row']), q, 1e-6)
            ok &= len(r['exposure_indices']) == 17
            if model.endswith('1'):
                expected_shift = f"{r['order_offset_m_s']:.2f}\\pm{r['order_offset_sigma_m_s']:.2f}"
                ok &= expected_shift in source and f'{c0["chi2"] - q:.3f}' in source
            if not ok:
                failure.append(f'{line}/{model}')
            rows_ok &= ok
    record('Order table reproduces eight matched native-pixel fits and additive Q', rows_ok, failure)
    later = s['2600']['heldout2019_2020']
    # Saved transfer objective is conditional; no independence or p-value inferred.
    dg = later['G0']['chi2'] - later['G1']['chi2']
    da = later['A1']['chi2'] - later['G1']['chi2']
    record('Conditional transfer improvements match prose', f'{dg:.2f}' in source and f'{da:.2f}' in source,
           {'G0_minus_G1': dg, 'A1_minus_G1': da})
    support = read(EXP / 'results/order_response/residual_common_support_summary.json')
    points = support['points']
    values = {}
    for line in [2382, 2600]:
        selected = [p for p in points if p['line'] == line]
        w = np.asarray([1 / p['variance_m2_s2'] for p in selected])
        d = np.asarray([p['delta_m_s'] for p in selected])
        values[str(line)] = [float(w @ d / w.sum()), float(1 / np.sqrt(w.sum()))]
    ok = all(allclose(values[k], [support['summary'][k]['mean_m_s'], support['summary'][k]['sigma_m_s']])
             for k in values)
    ok &= support['total_retained_good'] == sum(r['retained'] for r in support['retained_by_line'].values())
    ok &= support['total_original_good'] == sum(r['original'] for r in support['retained_by_line'].values())
    record('Common-support means, errors and native-pixel accounting recompute', ok, values)


def verify_cross_age():
    s = read(EXP / 'results/cross_age/design.json')
    z = np.array([r['z_abs'] for r in s['design_targets']])
    age = np.array([r['cosmic_age_Gyr'] for r in s['design_targets']])
    g = np.array([r['G2_z3_normalized'] for r in s['design_targets']])
    bg = s['background']
    h0_s = bg['H0_km_s_Mpc'] / 3.0856775814913673e19
    om = bg['Omega_m']
    def age_s(z):
        return 2 / (3*h0_s*np.sqrt(1-om)) * np.arcsinh(np.sqrt((1-om)/om) / (1+z)**1.5)
    sec_gyr = 365.25 * 86400 * 1e9
    t0, t3 = age_s(0), age_s(3)
    def G(p):
        return ((np.log(t0/bg['tstar_s']) / np.log(age_s(z)/bg['tstar_s']))**p - 1) / ((np.log(t0/bg['tstar_s']) / np.log(t3/bg['tstar_s']))**p - 1)
    record('Six age coordinates and normalized inverse-log response recompute independently',
           len(z) == 6 and allclose(age, age_s(z)/sec_gyr, 1e-6) and allclose(g, G(2), 1e-9))
    expected_shapes = {}
    for key, response in [('inverse_log_p1', G(1)), ('inverse_log_p3', G(3))]:
        design = np.column_stack([np.ones(len(z)), response])
        residual = g - design @ np.linalg.lstsq(design, g, rcond=None)[0]
        expected_shapes[key] = float(100*np.sqrt(np.mean(residual**2)))
    record('Conditional p1/p2/p3 shape-distance arithmetic recomputes',
           all(close(v, s['competitor_shapes'][k]['rms_residual_for_hypothetical_A_100_m_s'], 1e-8)
               for k, v in expected_shapes.items()), expected_shapes)
    delta = np.array([r['observed_pair_separation_AA'] for r in s['design_targets']])
    design = np.column_stack([np.ones(len(z)), delta/1000])
    beta = np.linalg.lstsq(design, 100*g, rcond=None)[0]
    residual = 100*g - design @ beta
    mimic = s['shared_slope_mimic_of_hypothetical_A_100_m_s']
    record('Shared-slope hypothetical magnitude recomputes',
           close(beta[1], mimic['fitted_shared_velocity_slope_m_s_per_1000_AA'], 1e-6)
           and close(np.sqrt(np.mean(residual**2)), mimic['residual_rms_m_s']))
    record('Independent differential nuisances span response, with no observed amplitude fit',
           np.linalg.matrix_rank(np.diag(delta)) == 6
           and s['no_data_dependent_amplitude_selection']
           and s['nuisance_projection']['one_wavelength_slope_per_target']['structural_nonidentifiability'])
    bad = hash_mapping(s['source_hashes'], EXP)
    record('Cross-age fixed-design source checksums match', not bad, bad)


def verify_archives():
    # The final six-target schema is checked here after final export exists.
    p = EXP / 'results/archive_complete/six_target_comparison.json'
    s = read(p)
    rows = s.get('comparisons', s.get('targets', s.get('rows')))
    if isinstance(rows, dict):
        rows = list(rows.values())
    if rows is None:
        raise KeyError('No comparisons/targets/rows in six-target summary')
    expected = {'J232128-105122', 'J225719-100104', 'J053007-250329',
                'J233156-090802', 'J004131-493611', 'J064326-504112'}
    record('All six screened archive targets remain in final comparison',
           len(rows) == 6 and {r['target'] for r in rows} == expected)
    table = text(PAPER / 'archive_table.tex')
    record('All six archive targets remain in manuscript table',
           all(t in table.replace('$-$', '-').replace('−', '-') for t in expected))
    for r in rows:
        n = read(EXP / r['null_fit_file'])
        a = read(EXP / r['alternative_fit_file'])
        k0 = n.get('npar', n.get('npar_total'))
        k1 = a.get('npar', a.get('npar_total'))
        matched = (n['ndata'] == a['ndata'] == r['pixels']
                   and n['ncomp'] == a['ncomp'] == r['ncomp']
                   and k1-k0 == r['extra_shifts']
                   and r['pixels']-k0 == r['nu0']
                   and r['pixels']-k1 == r['nu1']
                   and close(n['chi2'], r['Q0']) and close(a['chi2'], r['Q1'])
                   and close(n['chi2']-a['chi2'], r['delta_Q']))
        record(r['target'] + ': final archive Q, dimension and matched-fit accounting', matched)
        source_identity = sha(EXP / r['source_summary']) == r['source_summary_sha256']
        gate = (r['Q0']/r['nu0'] <= 1.5
                and max(x['chi2_per_pixel'] for x in r['null_per_line']) <= 1.8
                and n['optimizer_success'])
        d = r['D_m_s'] if r['adequacy_flag'] else r['diagnostic_only_D_m_s']
        source_identity &= close(d, a['shifts_m_s']['2382']) and gate == r['adequacy_flag']
        source_identity &= not r['physical_precision_measurement']
        source_identity &= not r['first_order_stationarity_certified']
        record(r['target'] + ': source identity, descriptor and adequacy bookkeeping', source_identity)
        target = r['target'].replace('-', '$-$')
        status = 'Conditional' if r['adequacy_flag'] else 'Inadequate'
        fragment = (f"{target} ({r['z_abs']:.3f}) & {r['pixels']} & {r['ncomp']} & "
                    f"{r['Q0']:.2f} & {r['delta_Q']:.2f} & {r['extra_shifts']} & {status}")
        record(r['target'] + ': rounded manuscript table row matches', fragment in table)
    stationarity = read(EXP / 'results/archive_complete/validation/stationarity_independent.json')
    binding_errors = []
    for r in rows:
        if 'common_pair/' in r['null_fit_file']:
            continue
        for key in ['null_fit_file', 'alternative_fit_file']:
            path = EXP / r[key]
            fit = read(path)
            matches = [d for d in stationarity['fit_diagnostics']
                       if d['target'] == r['target'] and d['name'] == fit['name']]
            if len(matches) != 1:
                binding_errors.append(r[key])
                continue
            d = matches[0]
            if (d['source_json_sha256'] != sha(path)
                    or d['source_npz_sha256'] != sha(path.with_suffix('.npz'))):
                binding_errors.append(r[key])
    record('Eight final new-archive fit artifacts match independent audit hashes',
           not binding_errors, binding_errors)
    r225 = next(r for r in rows if r['target'] == 'J225719-100104')
    selection = read(EXP / 'results/archive_complete/validation/per_count_selection_audit.json')
    gaps = selection['lower_incomplete_gaps']
    record('Lower unfinished J2257 endpoint and recovered comparison remain visible',
           len(gaps) >= 1 and any(r['target'] == 'J225719-100104' for r in gaps)
           and r225['Q0'] < min(r['lowest_incomplete_chi2'] for r in gaps if r['target'] == 'J225719-100104')
           and any(not r['optimizer_success'] for r in r225['sparse_profile']))
    record('Archive export source identity recorded by table provenance',
           read(PAPER / 'table_provenance.json')['experiment_source_sha256'].get(
               'results/archive_complete/six_target_comparison.json') == sha(p))


def verify_validation_reports():
    cases = [
        ('results/common_pair/review/validation.json', 'passed'),
        ('results/common_pair/review/delivery_validation.json', 'passed'),
        ('results/order_response/validation.json', 'pass_'),
        ('results/order_response/residual_validation.json', 'pass'),
        ('results/cross_age/independent_validation.json', 'status'),
        ('results/archive_complete/validation/independent_validation.json', 'passed'),
        ('results/archive_complete/validation/stationarity_independent.json', 'all_checks_pass'),
    ]
    for rel, key in cases:
        s = read(EXP / rel)
        value = s[key]
        ok = value == 'PASS' if key == 'status' else value is True
        ok &= not s.get('partial', False)
        record('Completed source audit: ' + rel, ok,
               {'scope': s.get('scope'), 'checks': s.get('n_checks', s.get('checks_total', s.get('check_count')))})
    p = EXP / 'data/raw/calibration_readiness/source_hash_manifest.json'
    s = read(p)
    bad = []
    for row in s['files']:
        path = EXP / row['path']
        if not path.is_file() or sha(path) != row['sha256'] or path.stat().st_size != row['bytes']:
            bad.append(row['path'])
    record('Calibration-readiness artifact hash and byte manifest matches',
           not bad and len(s['files']) == s['file_count'], bad)
    c = read(EXP / 'data/raw/calibration_readiness/validation.json')
    record('Calibration-readiness scoped audit has no failed checks',
           c['passed'] == c['total'] and all(r['pass'] for r in c['checks']))
    lit = read(ROOT / 'data/raw/spectroscopy_literature/manifest.json')
    record('New primary-literature metadata checksums match',
           all(sha(ROOT / f'data/raw/spectroscopy_literature/{r["key"]}_crossref.json') == r['sha256'] for r in lit))


def verify_markdown():
    paths = [ROOT / 'readme.md', ROOT / 'readme_cn.md', EXP / 'readme.md',
             EXP / 'reports/publication_followup_results_cn.md']
    missing = []
    links = 0
    for path in paths:
        content = text(path)
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)', content):
            target = target.strip().strip('<>')
            if '://' in target or target.startswith(('mailto:', '#')):
                continue
            target = unquote(target.split('#')[0])
            # Accept line suffix only when it is a local file-link annotation.
            target = re.sub(r':\d+$', '', target)
            if not target:
                continue
            links += 1
            p = Path(target) if target.startswith('/') else path.parent / target
            if not p.exists():
                missing.append({'from': str(path.relative_to(ROOT)), 'target': target})
    record('Local links in both project READMEs, experiment README and new summary resolve',
           not missing, {'links_checked': links, 'missing': missing})
    record('Both project READMEs link the focused manuscript',
           all('paper_spectroscopy/feii_relative_frequencies.pdf' in (ROOT / n).read_text()
               for n in ['readme.md', 'readme_cn.md']))


def verify_pdf():
    pdf = PAPER / 'feii_relative_frequencies.pdf'
    p = read(REPORT / 'spectroscopy_pdf_preflight.json')
    record('Final PDF preflight passes and matches delivered PDF hash',
           p['verdict'] == 'PASS' and p['sha256'] == sha(pdf)
           and p['enumerated_page_count'] == p['declared_page_count'] == p['reader_page_count'],
           {'pages': p.get('reader_page_count'), 'sha256': sha(pdf)})
    record('Delivered PDF matches current build PDF', sha(pdf) == sha(BUILD / 'main.pdf'))
    log = text(BUILD / 'main.log')
    bad = re.findall(r'[^\n]*(?:undefined|multiply defined|LaTeX Error|Emergency stop|Fatal error)[^\n]*', log, flags=re.I)
    record('Final LaTeX log has no undefined references or fatal errors', not bad, bad)
    biblog = text(BUILD / 'main.blg')
    errors = re.findall(r'[^\n]*(?:Warning--I didn.t find|I couldn.t open|error message)[^\n]*', biblog)
    record('BibTeX found cited entries and bibliography input', not errors, errors)
    newest = max(p.stat().st_mtime for p in PAPER.glob('*.tex'))
    newest = max(newest, (PAPER / 'references.bib').stat().st_mtime)
    record('PDF was built after latest manuscript source modification', pdf.stat().st_mtime >= newest)


def verify_prediction_appendix():
    prediction = ROOT / 'experiments/prediction_test_2026-09-21'
    s = read(prediction / 'results/design.json')
    appendix = text(PAPER / 'prediction_appendix.tex')
    body = text(PAPER / 'discussion.tex')
    main_source = text(PAPER / 'main.tex')
    bg = s['cosmology']
    z = np.array([r['z'] for r in s['geometry']])
    age = 2 / (3*np.sqrt(1-bg['Omega_m'])) * np.arcsinh(
        np.sqrt((1-bg['Omega_m'])/bg['Omega_m'])/(1+z)**1.5)
    age_seconds = age * 3.0856775814913673e19 / bg['H0_km_s_Mpc']
    raw = 1 / np.log(age_seconds/bg['tstar_seconds'])**2
    fraction = (raw[1]-raw[0])/(raw[-1]-raw[0])
    weights = np.array([1-fraction, fraction])
    stored = s['conditional_prediction']
    record('Prospective interpolation independently matches fixed inverse-log-square law',
           allclose(weights, stored['anchor_weights'], 1e-10)
           and close(100*fraction, stored['example_middle_minus_low_m_s'], 1e-8))
    record('Brief prospective section and detailed appendix are included',
           r'\label{sec:prospective}' in body and r'\label{app:prediction}' in appendix
           and r'\input{prediction_appendix.tex}' in main_source
           and r'\appendix' in main_source)
    record('Main-text prospective coefficients and example match design rounding',
           f'{weights[0]:.5f}' in body and f'{weights[1]:.5f}' in body
           and f'{100*fraction:.2f}' in body)
    record('Prospective design remains hypothetical and internally checked',
           stored['example_is_not_measured_or_theory_predicted']
           and s['validation']['passed']
           and not any(s['operational_status'].values())
           and all(c['passed'] for c in s['validation']['checks']))
    record('Prospective code and atomic source identities match frozen design',
           s['source_hashes']['design_code'] == sha(prediction/'code/prediction_design.py')
           and s['source_hashes']['atomic_identification_csv'] == sha(EXP/'data/atomic/selected_transitions.csv'))
    figure = prediction/'figures/prediction_design.pdf'
    INPUTS[str(figure.relative_to(ROOT))] = sha(figure)
    read(prediction/'reports/source_manifest.json')
    text(prediction/'reports/instrument_feasibility_cn.md')
    text(prediction/'reports/prediction_logic_review.md')
    text(prediction/'prediction_protocol_cn.md')


def main():
    for name, fn in [('TeX', verify_tex), ('Control export', verify_controls),
                     ('Common pair', verify_common_pair), ('Order response', verify_order),
                     ('Cross-age design', verify_cross_age), ('Archive final export', verify_archives),
                     ('Source audits and manifests', verify_validation_reports),
                     ('README and summary links', verify_markdown),
                     ('Prospective appendix', verify_prediction_appendix), ('PDF delivery', verify_pdf)]:
        audit(name, fn)
    failed = [r for r in CHECKS if not r['passed']]
    out = {
        'scope': 'Focused-manuscript delivery integrity, selected number recomputation, exported table/source identity, local link and final build audit only. No new nonlinear fits, significance calibration, coverage guarantee, physical detection certification or publication acceptance claim.',
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'passed': not failed,
        'checks_passed': len(CHECKS)-len(failed),
        'checks_total': len(CHECKS),
        'checks': CHECKS,
        'failures': failed,
        'input_sha256': INPUTS,
        'validator_sha256': sha(Path(__file__)),
    }
    REPORT.mkdir(exist_ok=True)
    (REPORT / 'spectroscopy_delivery_verification.json').write_text(json.dumps(out, indent=2, ensure_ascii=False)+'\n')
    lines = ['# Focused spectroscopy delivery audit', '',
             f"Status: **{'PASS' if out['passed'] else 'INCOMPLETE / FAIL'}** ({out['checks_passed']}/{out['checks_total']} grouped checks).", '',
             out['scope'], '',
             'These grouped delivery checks do not add to the individual scientific validation counts reported by the source analyses.', '',
             '| Check | Result |', '|---|---|']
    lines += [f"| {r['name']} | {'PASS' if r['passed'] else 'FAIL'} |" for r in CHECKS]
    if failed:
        lines += ['', 'Unresolved items:', '']
        lines += [f"- {r['name']}: `{str(r.get('detail',''))}`" for r in failed]
    lines += ['', 'The JSON companion records exact input checksums and details.']
    (REPORT / 'spectroscopy_delivery_verification.md').write_text('\n'.join(lines)+'\n')
    print(f"{'PASS' if out['passed'] else 'INCOMPLETE / FAIL'}: {out['checks_passed']}/{out['checks_total']} grouped checks")
    for r in failed:
        print(r['name'], r.get('detail', ''))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
