#!/usr/bin/env python3
"""Continue the recorded capped core01 branch and compare matched endpoints.

The original attempt is retained. Existing continuation artifacts are reused;
delete only this script's named products to repeat the optimization itself.
Run after espresso_saturation_controls.py has written both comparisons.
"""
import json
from espresso_saturation_controls import fit, OUT


def get_or_fit(name, free, start):
    path = OUT / (name + '.json')
    if path.exists():
        return json.loads(path.read_text())
    return fit(name, free, .1, 0., OUT / (start + '.json'), 500)


def main():
    a = get_or_fit('core01_alternative_cross_continued', True, 'core01_alternative_cross')
    get_or_fit('core01_null_cross_continued_alternative', False, a['name'])
    paths = [p for p in OUT.glob('core01_*.json') if not any(t in p.name for t in ['checkpoint', 'comparison'])]
    trials = [json.loads(p.read_text()) for p in paths]
    assert all(r['ndata'] == 2692 for r in trials)
    # Preserve unsuccessful endpoints too, so their objective values are visible.
    null = min((r for r in trials if not r['free_shifts'] and r['optimizer_success']), key=lambda r:r['chi2'])
    alt = min((r for r in trials if r['free_shifts'] and r['optimizer_success']), key=lambda r:r['chi2'])
    unresolved = [r['name'] for r in trials if not r['optimizer_success'] and
                  r['chi2'] < (alt if r['free_shifts'] else null)['chi2']]
    assert not unresolved, f'Lower unfinished endpoints must be disclosed or continued: {unresolved}'
    result = dict(case='core01', null=null['name'], alternative=alt['name'], ndata=null['ndata'],
                  added_parameters=alt['npar']-null['npar'], delta_chi2=null['chi2']-alt['chi2'],
                  null_success=null['optimizer_success'], alternative_success=alt['optimizer_success'],
                  shifts_m_s=alt['shifts_m_s'], mask_audit=null['mask_audit'],
                  selection='Lowest objective among successful recorded endpoints; optimizer termination does not certify stationarity or global minimality.',
                  recorded_trials=[{k:r[k] for k in ['name','free_shifts','chi2','optimizer_success','nfev','optimality']} for r in trials])
    (OUT / 'core01_comparison.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
