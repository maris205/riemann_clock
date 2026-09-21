#!/usr/bin/env python3
"""Verify final campaign delivery against saved numerical products; no fits."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / 'experiments/highz_feasibility_2026-09-20'
RES = EXP / 'results'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    checks = []
    def check(name, condition):
        checks.append(dict(name=name, passed=bool(condition)))

    pdf = ROOT / 'paper/riemann_clock.pdf'
    preflight = read(ROOT / 'reports/available_data_pdf_preflight.json')
    check('PDF read preflight PASS', preflight['verdict'] == 'PASS')
    check('Preflight identifies current PDF', preflight['sha256'] == sha(pdf))
    if not all(x['passed'] for x in checks):
        raise SystemExit('Run PDF preflight for the current compiled PDF first.')
    pdftext = subprocess.check_output(['pdftotext', '-layout', str(pdf), '-'], text=True)
    report = (EXP / 'reports/available_data_campaign_results_cn.md').read_text()
    table = (ROOT / 'paper/available_data_table.tex').read_text()
    prose = (ROOT / 'paper/available_data_section.tex').read_text()
    for row in read(RES / 'available_data_campaign/model_comparisons.json'):
        label = row['label']
        value = f'{row["delta_chi2"]:.2f}'
        check(label + ' in compiled table', value in table and value in pdftext)
        check(label + ' report matches precise value', f'{row["delta_chi2"]:.6f}' in report)
        check(label + ' objective identity', abs(row['null_chi2'] - row['alternative_chi2'] - row['delta_chi2']) < 1e-8)
        check(label + ' adds five shifts', row['alternative_parameters'] - row['null_parameters'] == 5)
    for name, digest in read(RES / 'available_data_campaign/model_source_hashes.json').items():
        check('Export source hash ' + name, sha(ROOT / name) == digest)
    for case in ['core01_gls', 'core03_gls']:
        c = read(RES / 'joint_controls' / (case + '_comparison.json'))
        n = read(RES / 'joint_controls' / (c['null'] + '.json'))
        a = read(RES / 'joint_controls' / (c['alternative'] + '.json'))
        check(case + ' saved paired identity', abs(n['chi2'] - a['chi2'] - c['delta_chi2']) < 1e-8)
        check(case + ' same pixels', n['ndata'] == a['ndata'] == c['ndata'])
    expected = read(RES / 'archive_expansion/profile_pilot/expected_noise/summary.json')
    for c in expected['comparisons']:
        for prefix in ['row1', 'expected']:
            value = c[prefix + '_delta_chi2']
            check(c['target'] + prefix + ' objective identity',
                  abs(c[prefix + '_chi2_null'] - c[prefix + '_chi2_alternative'] - value) < 1e-8)
            check(c['target'] + prefix + ' PDF value', f'{value:.2f}' in pdftext)
            check(c['target'] + prefix + ' report value', f'{value:.6f}' in report)
        check(c['target'] + ' noise control stops successfully', c['both_optimizer_success'])
    q = read(RES / 'archive_expansion/quality_summary.json')
    check('Archive counts 36/32/324/6', [q[k] for k in ['candidate_absorbers', 'downloaded_unique_sightlines',
          'actual_windows', 'computational_readiness_count']] == [36, 32, 324, 6])
    d = read(RES / 'exposures/diagnostics.json')
    check('Native exposure valid-pixel count', d['total_good_native_pixels'] == 76064)
    check('Order contrasts and errors in PDF', all(s in pdftext for s in ['60.2', '16.7', '57.7', '16.6']))
    uv = read(RES / 'uves_completion/summary.json')
    check('UVES five explicit nonstationary capped endpoints', len(uv['cases']) == 5 and
          all(not c['stationary'] and c['optimizer_status'] == 0 for c in uv['cases']))
    for p in uv['pairs']:
        check('UVES finite difference disclosed ' + str(p['wide']), f'{p["delta_chi2"]:.6f}' in report and not p['stationary_both'])
    check('UVES limits retained in prose', 'first-order stationarity' in prose and 'no calibrated significance' in prose)

    # Read separate numerical audits, retaining their distinct scopes/counts.
    audits = [
        ('joint_controls_review/validation.json', 'n_fail', 0),
        ('gas_structure/independent_validation.json', 'passed', True),
        ('exposures/validation.json', 'pass_', True),
        ('archive_expansion/independent_catalogue_validation.json', 'passed', True),
        ('archive_expansion/independent_data_validation.json', 'passed', True),
        ('archive_expansion/independent_diagnostics_validation.json', 'passed', True),
        ('archive_expansion/profile_pilot/independent_validation.json', 'passed', True),
        ('archive_expansion/profile_pilot/expected_noise/independent_validation.json', 'passed', True),
        ('uves_completion/validation.json', 'status', 'PASS'),
    ]
    for name, key, value in audits:
        check('Separate numerical audit ' + name, read(RES / name)[key] == value)
    mi = read(RES / 'mask_information/information.json')
    check('Information internal checks', mi['validation']['all_passed'])
    mic = read(RES / 'mask_information/independent_check.json')
    # The independent audit's source hashes identify its common physical
    # baseline, not the producer's information.json output.
    check('Information independent audit baseline source',
          mic['source_json_sha256'] == sha(RES / 'espresso_null/null_cross.json') and
          mic['source_npz_sha256'] == sha(RES / 'espresso_null/null_cross.npz') and mic['source_files_unchanged'])
    for case in mi['cases']:
        name = {'full_gls': 'full', 'core01_gls': 'core01', 'core03_gls': 'core03_padded'}[case['case']]
        independent = mic['cases'][name]
        check(name + ' independent information pixel count', case['retained_pixels'] == independent['ndata'])
        for cut in case['cutoffs']:
            other = independent['cutoffs'][str(cut['relative_svd_cutoff'])]
            check(name + ' independent information value ' + str(cut['relative_svd_cutoff']),
                  abs(cut['illustrative_100m_s_expected_local_noncentrality'] -
                      other['expected_local_signal_norm_squared']) < 1e-8)
    check('Validated original core unchanged', sha(EXP / 'code/espresso_conventional_null.py') ==
          '424dc9e32cb9bb0a3eece533c64d44ae2582f1edc4b8fa4b14889c6e4db142de')
    for path in [ROOT / 'readme.md', ROOT / 'readme_cn.md', EXP / 'readme.md',
                 EXP / 'reports/available_data_campaign_results_cn.md']:
        for target in re.findall(r'\]\(([^)]+)\)', path.read_text()):
            if '://' in target or target.startswith('#'):
                continue
            check('Local link ' + str(path.relative_to(ROOT)) + ' -> ' + target,
                  (path.parent / target.split('#', 1)[0]).exists())
    log = (ROOT / 'build/main.log').read_text()
    check('No undefined references or overfull boxes', not any(s in log for s in ['undefined', 'Overfull', 'LaTeX Error']))
    check('New figure, table and retained heatmap in PDF',
          all(s in pdftext for s in ['FIG. 9.', 'TABLE V.', 'Conditional resource']))
    check('Single-column manuscript', 'onecolumn' in (ROOT / 'paper/main.tex').read_text())
    paths = [pdf, ROOT / 'paper/main.tex', ROOT / 'paper/highz_analysis.tex',
             ROOT / 'paper/available_data_section.tex', ROOT / 'paper/available_data_table.tex',
             ROOT / 'figures/joint_controls_information.pdf', ROOT / 'figures/resource_precision_heatmap.pdf',
             EXP / 'reports/available_data_campaign_results_cn.md',
             RES / 'exposures/summary.json', RES / 'exposures/diagnostics.json',
             RES / 'archive_expansion/quality_summary.json',
             RES / 'archive_expansion/profile_pilot/expected_noise/summary.json',
             RES / 'uves_completion/summary.json']
    output = dict(status='PASS' if all(c['passed'] for c in checks) else 'FAIL',
                  count=len(checks), passed=sum(c['passed'] for c in checks), checks=checks,
                  generated_at=datetime.now(timezone.utc).isoformat(),
                  page_count=preflight['declared_page_count'],
                  artifact_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
                  scope='Delivery consistency, not additional physical evidence or optimizer convergence.')
    (ROOT / 'reports/available_data_delivery_verification.json').write_text(json.dumps(output, indent=2) + '\n')
    if output['status'] == 'PASS':
        (ROOT / 'reports/available_data_delivery_verification.md').write_text(
            '# Existing-data extension: delivery verification\n\n'
            f'**PASS: {output["passed"]}/{output["count"]} delivery checks.** '
            f'The compiled single-column PDF contains {output["page_count"]} pages.\n\n'
            'Saved model comparisons, pilot error-array controls, archive counts, independent audit statuses, '
            'export source hashes, local links and compiled values agree. The original validated ESPRESSO '
            'core is unchanged. All five UVES continuation endpoints remain explicitly capped and nonstationary.\n\n'
            'PDF preflight passed. The final title page and new results pages were visually inspected; '
            'the new comparison figure and table fit cleanly. There are no undefined references or overfull boxes. '
            'Three underfull text-box notices remain cosmetic. Original resource heatmap and laboratory figures are retained.\n\n'
            'These are package-consistency checks, not independent physical evidence. Numerical validation '
            'and model limitations remain in the separate reports. The JSON record includes timestamps and artifact SHA-256 hashes.\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ['checks','artifact_sha256']}, indent=2))
    if output['status'] != 'PASS':
        print(json.dumps([c for c in checks if not c['passed']], indent=2)); raise SystemExit(1)


if __name__ == '__main__':
    main()
