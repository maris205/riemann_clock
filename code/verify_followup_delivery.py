#!/usr/bin/env python3
"""Check the final follow-up package, without repeating numerical fits."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / 'experiments/highz_feasibility_2026-09-20'
checks = []


def check(name, condition):
    checks.append(dict(name=name, passed=bool(condition)))


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    summary = read(EXP / 'results/followup_summary/summary.json')
    table = (ROOT / 'paper/espresso_followup_table.tex').read_text()
    macros = (ROOT / 'paper/espresso_followup_numbers.tex').read_text()
    report = (EXP / 'reports/conventional_followup_results_cn.md').read_text()
    pdftext = (ROOT / 'build/followup_pdf_text.txt').read_text()
    for row in summary['espresso_comparisons']:
        value = f'{row["delta_chi2"]:.2f}'
        check(row['case'] + ' number in TeX table', value in table)
        check(row['case'] + ' number in Chinese report', value in report)
        check(row['case'] + ' number in compiled PDF', value in pdftext)
        if row['case'] != 'primary':
            check(row['case'] + ' macro', '{' + value + '}' in macros)
            a = read(EXP / row['null_file']); b = read(EXP / row['alternative_file'])
            check(row['case'] + ' final endpoint identity', abs(a['chi2'] - b['chi2'] - row['delta_chi2']) < 1e-9)
            check(row['case'] + ' paired pixels and parameter count', a['ndata'] == b['ndata'] and b['npar'] - a['npar'] == 5)
    for name, digest in read(EXP / 'results/followup_summary/source_hashes.json').items():
        check('source hash ' + name, sha(ROOT / name) == digest)
    preflight = read(ROOT / 'reports/followup_manuscript_pdf_preflight.json')
    check('PDF preflight', preflight['verdict'] == 'PASS' and preflight['declared_page_count'] == 32)
    check('PDF preflight refers to final PDF', preflight['sha256'] == sha(ROOT / 'paper/riemann_clock.pdf'))
    log = (ROOT / 'build/main.log').read_text()
    check('no undefined references or overfull boxes', not any(t in log for t in ['Overfull', 'undefined', 'LaTeX Error']))
    check('resource heatmap retained', 'resource_precision_heatmap.pdf' in (ROOT / 'paper/main.tex').read_text() and 'Conditional resource' in pdftext)
    check('new figure included', 'FIG. 8.' in pdftext)
    check('HARPERFECT citation compiled', 'HARPS spectrograph' in pdftext)
    check('UVES limited optimization disclosed', 'Both null searches' in pdftext and '500 evaluations' in pdftext)
    for path in [ROOT/'readme.md', ROOT/'readme_cn.md', EXP/'readme.md', EXP/'reports/conventional_followup_results_cn.md']:
        for target in re.findall(r'\]\(([^)]+)\)', path.read_text()):
            if '://' in target or target.startswith('#'):
                continue
            clean = target.split('#', 1)[0]
            check('local link ' + str(path.relative_to(ROOT)) + ' -> ' + target, (path.parent / clean).exists())
    for path in [EXP/'results/noise_covariance/independent_validation.json', EXP/'results/cross_instrument/validation.json']:
        check(str(path.relative_to(ROOT)) + ' PASS', read(path)['status'] == 'PASS')
    reviewed = read(EXP/'results/followup_review/validation.json')
    check('nonlinear independent review contains 512 checks', '512' in (EXP/'reports/espresso_followup_independent_review.md').read_text())
    check('nonlinear independent review all checks pass', all(c['passed'] for c in reviewed['checks']))
    artifact_paths = [ROOT/'paper/riemann_clock.pdf', ROOT/'paper/espresso_followup_section.tex',
                      ROOT/'figures/espresso_noise_controls.pdf', ROOT/'figures/resource_precision_heatmap.pdf',
                      EXP/'reports/conventional_followup_results_cn.md', EXP/'reports/cross_instrument_results_cn.md',
                      EXP/'results/cross_instrument/summary.json', EXP/'results/followup_summary/summary.json']
    result = dict(status='PASS' if all(c['passed'] for c in checks) else 'FAIL',
                  checks_passed=sum(c['passed'] for c in checks), count=len(checks), checks=checks,
                  generated_at=datetime.now(timezone.utc).isoformat(),
                  artifact_sha256={str(p.relative_to(ROOT)):sha(p) for p in artifact_paths},
                  scope='Delivery consistency only; numerical audits and explicitly unresolved UVES optimization remain separate.')
    (ROOT/'reports/followup_delivery_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['checks','artifact_sha256']}, indent=2))
    if result['status'] != 'PASS':
        print(json.dumps([c for c in checks if not c['passed']], indent=2))
        raise SystemExit(1)


if __name__ == '__main__':
    main()
