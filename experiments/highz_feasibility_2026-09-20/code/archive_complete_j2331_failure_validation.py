#!/usr/bin/env python3
"""Independent, no-fit replay of the bounded J2331 failure-case diagnostic."""
from pathlib import Path
import hashlib
import json
import numpy as np
from archive_complete_j2331_model import NeighborModel

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / 'results/archive_complete/J233156-090802'
OUT = PRIMARY / 'failure_pair'
CHECKS = []


def check(name, condition, detail=''):
    CHECKS.append((name, bool(condition), detail))


def close(name, actual, expected, rtol=2e-12, atol=2e-10):
    a, b = np.asarray(actual), np.asarray(expected)
    good = a.shape == b.shape and np.allclose(a, b, rtol=rtol, atol=atol)
    detail = 'shape mismatch' if a.shape != b.shape else f'max abs difference {np.max(np.abs(a-b)) if a.size else 0:.6g}'
    check(name, good, detail)


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate():
    protocol = read(OUT / 'protocol.json')
    summary = read(OUT / 'summary.json')
    cfg = read(PRIMARY / 'selection.json')
    original_path = ROOT / protocol['initial_primary_null_file']
    original = read(original_path)
    fitpaths = sorted((OUT / 'fits/J233156-090802').glob('*.json'))
    fitpaths = [p for p in fitpaths if not p.stem.endswith('_checkpoint')]
    paths = [original_path] + fitpaths
    records = {p.stem: read(p) for p in paths}
    arrays = {p.stem: np.load(p.with_suffix('.npz')) for p in paths}
    source = np.load(ROOT / 'data/processed/archive_expansion/J233156-090802.npz')
    checks_by_record = {}

    check('protocol config frozen hash', sha(PRIMARY / 'selection.json') == protocol['configuration_sha256'])
    check('summary config hash', sha(PRIMARY / 'selection.json') == summary['primary_configuration_sha256'])
    check('summary protocol hash', sha(OUT / 'protocol.json') == summary['protocol_sha256'])
    check('protocol configuration contents', protocol['configuration'] == cfg)
    check('original primary JSON unchanged', sha(original_path) == protocol['initial_primary_null_sha256'])
    check('original primary NPZ unchanged', sha(original_path.with_suffix('.npz')) == protocol['initial_primary_null_npz_sha256'])
    for field, path in [
        ('model_adapter_sha256', ROOT / 'code/archive_complete_j2331_model.py'),
        ('shared_engine_sha256', ROOT / 'code/archive_complete_model.py'),
        ('runner_sha256', ROOT / 'code/archive_complete_j2331_failure_pair.py')]:
        check('protocol ' + field, sha(path) == protocol[field])
    close('recorded original objective', original['chi2'], protocol['initial_primary_null_chi2'])
    check('profile uncertainty disabled', protocol['profile_uncertainty'] is False)
    check('fixed shift bounds protocol', protocol['free_shift_bounds_km_s'] == [-1., 1.])
    check('original precision disposition preserved', summary['primary_precision_descriptor_D'] == 'not estimated; unchanged')
    check('diagnostic marked nonprecision', summary['diagnostic_only'] is True and summary['primary_adequacy_gate'] is False)
    ordinary = read(PRIMARY / 'ordinary_model_summary.json')
    check('original gate still false', ordinary['passes_conditional_adequacy_gate'] is False)

    common_ids = None
    models = {}
    hashes = {}
    sourcefiles = {
        'coadd': ROOT / 'data/processed/archive_expansion/J233156-090802.npz',
        'metadata': ROOT / 'data/processed/archive_expansion/J233156-090802_metadata.json',
        'atomic': ROOT / 'data/raw/espresso_null/MM_VPFIT_2013-11-10.dat',
        'base_model': ROOT / 'code/archive_profile_pilot.py',
        'adapter': ROOT / 'code/archive_complete_model.py'}
    for path in paths:
        before = len(CHECKS)
        r, a = records[path.stem], arrays[path.stem]
        prefix = path.stem + ': '
        model = NeighborModel(cfg, 26, r['free_shifts'], 21)
        models[path.stem] = model
        p = np.array(r['parameters'])
        _, lo, hi = model.initial(r['seed'])
        res, jac, profiles = model.evaluate(p)
        check(prefix + 'fixed architecture', r['ncomp'] == 26 and r['oversample'] == 21 and r['ndata'] == 1176)
        check(prefix + 'same configuration', r['configuration'] == cfg)
        check(prefix + 'labels unchanged', r['labels'] == model.labels)
        check(prefix + 'parameter count', r['npar'] == model.npar == (92 if r['free_shifts'] else 90))
        check(prefix + 'no physical alpha or time parameter', not any('alpha' in s or 'time' in s for s in model.labels))
        check(prefix + 'bounds observed', np.all(p >= lo) and np.all(p <= hi))
        close(prefix + 'lower bounds', r['bounds_lower'], lo)
        close(prefix + 'upper bounds', r['bounds_upper'], hi)
        close(prefix + 'NPZ parameters', a['parameters'], p)
        close(prefix + 'residual replay', a['residual'], res)
        close(prefix + 'full Jacobian replay', a['jacobian'], jac)
        check(prefix + 'finite residual and Jacobian', np.isfinite(res).all() and np.isfinite(jac).all())
        check(prefix + 'all parameters fitted', a['fitted_parameter_mask'].shape == (model.npar,) and a['fitted_parameter_mask'].all())
        chi = float(res @ res)
        close(prefix + 'objective replay', r['chi2'], chi)
        check(prefix + 'nominal dof', r['nominal_ndf'] == 1176 - model.npar)
        close(prefix + 'reduced objective', r['chi2_per_ndf'], chi / r['nominal_ndf'])
        close(prefix + 'AICc arithmetic', r['aicc'], chi + 2*model.npar + 2*model.npar*(model.npar+1)/(1176-model.npar-1))
        check(prefix + 'termination status retained', r['optimizer_success'] == (r['optimizer_status'] > 0))
        check(prefix + 'evaluation budget respected', 0 < r['nfev'] <= r['max_nfev'])
        check(prefix + 'no fixed shift parameter', r['fixed_2382_shift_km_s'] is None)
        expected_active = [model.labels[i] for i in range(model.npar) if min(p[i]-lo[i], hi[i]-p[i]) < 1e-5]
        check(prefix + 'active bounds replay', expected_active == r['active_bounds'])
        for field, sourcepath in sourcefiles.items():
            check(prefix + 'source hash ' + field, sha(sourcepath) == r['source_hashes'][field])
        ids = []
        for L, sl, profile, savedline in zip(model.lines, model.slices, profiles, r['per_line']):
            key = str(L['key'])
            lp = prefix + key + ': '
            ix = L['source_indices']
            independent_v = 299792.458 * np.log(source['wavelength_AA']/(L['ref']*(1+cfg['z'])))
            independent_ix = np.flatnonzero((independent_v >= -480) & (independent_v <= 500))
            check(lp + 'independent window indices', np.array_equal(ix, independent_ix))
            check(lp + '392 valid samples', len(ix) == 392 and L['good'].all())
            for field in ['v', 'wave', 'flux', 'error', 'statistical_error', 'good', 'source_indices']:
                check(lp + field + ' exact saved equality', np.array_equal(a[key+'_'+field], L[field]))
            check(lp + 'raw flux', np.array_equal(L['flux'], source['flux'][ix]))
            check(lp + 'expected error row', np.array_equal(L['error'], source['expected_fluctuation'][ix]))
            close(lp + 'model flux replay', a[key+'_model'], profile)
            close(lp + 'direct saved-flux residual', (a[key+'_model'][L['good']]-source['flux'][ix][L['good']])/source['expected_fluctuation'][ix][L['good']], res[sl])
            check(lp + 'per-line count', savedline['line'] == L['key'] and savedline['ndata'] == 392)
            close(lp + 'per-line objective', savedline['chi2'], float(res[sl]@res[sl]))
            close(lp + 'per-line objective per pixel', savedline['chi2_per_pixel'], float(res[sl]@res[sl])/392)
            ids.extend(ix[L['good']].tolist())
        check(prefix + '1176 globally unique likelihood pixels', len(ids) == len(set(ids)) == 1176)
        if common_ids is None:
            common_ids = ids
        check(prefix + 'all records identical native pixels and order', ids == common_ids)
        check(prefix + 'no separate 1611 data region', not any(s.startswith('1611_') for s in a.files))
        g = jac.T @ res
        scale = hi-lo
        q = (p-lo)/scale
        cl = np.max(np.where(g > 0, (p-lo)*g, -(hi-p)*g))
        gm = np.max(np.abs(q-np.clip(q-g*scale, 0, 1)))
        cg = np.max(np.abs(g/np.maximum(np.linalg.norm(jac, axis=0), 1e-30)))
        close(prefix + 'independent Coleman Li', r['stationarity']['independent_coleman_li'], cl)
        close(prefix + 'normalized projected gradient', r['stationarity']['normalized_gradient_mapping_inf'], gm)
        close(prefix + 'column normalized gradient', r['stationarity']['column_normalized_gradient_inf'], cg)
        hashes[str(path.relative_to(ROOT))] = sha(path)
        hashes[str(path.with_suffix('.npz').relative_to(ROOT))] = sha(path.with_suffix('.npz'))
        checks_by_record[path.stem] = len(CHECKS)-before

    first = records['diagnostic_H1_from_primary_os21']
    cross = records['diagnostic_H0_cross_os21']
    need_cross = cross['chi2'] < original['chi2'] - .001
    extra = records.get('diagnostic_H1_cross_os21')
    check('optional cross trigger faithfully applied', (extra is not None) == need_cross)
    check('bounded number of new fits', len(fitpaths) == (3 if need_cross else 2))
    for name, start in [('diagnostic_H1_from_primary_os21', original), ('diagnostic_H0_cross_os21', first)] + ([('diagnostic_H1_cross_os21', cross)] if extra else []):
        r = records[name]
        model = models[name]
        p0, lo, hi = model.initial(r['seed'])
        mapped = dict(zip(start['labels'], start['parameters']))
        p0 = np.array([mapped.get(label, default) for label, default in zip(model.labels, p0)])
        close(name + ' exact initial mapping', r['initial_parameters'], np.clip(p0, lo+1e-8, hi-1e-8))
        check(name + ' saved after frozen protocol', (OUT/'protocol.json').stat().st_mtime <= next(p for p in paths if p.stem == name).stat().st_mtime)
    null = min([original, cross], key=lambda r: r['chi2'])
    alt = min([first]+([extra] if extra else []), key=lambda r: r['chi2'])
    check('selected retained H0', summary['null']['name'] == null['name'])
    check('selected retained H1', summary['alternative']['name'] == alt['name'])
    for label, rec in [('null', null), ('alternative', alt)]:
        check(label + ' compact summary copied exactly', all(value == rec[key] for key, value in summary[label].items()))
    close('delta objective arithmetic', summary['delta_chi2'], null['chi2']-alt['chi2'])
    check('two extra region shifts', summary['extra_region_shifts'] == alt['npar']-null['npar'] == 2)
    check('both stop flags preserved', summary['both_selected_optimizer_success'] == (null['optimizer_success'] and alt['optimizer_success']))
    gate = bool(alt['optimizer_success'] and alt['chi2_per_ndf'] <= 1.5 and max(r['chi2_per_pixel'] for r in alt['per_line']) <= 1.8)
    check('alternative numeric gate recomputed', summary['alternative_meets_same_numeric_quality_thresholds'] == gate)
    m1 = models[alt['name']]
    pembed = np.array([dict(zip(null['labels'], null['parameters'])).get(k, 0.) for k in m1.labels])
    er = m1.fun(pembed)
    close('H0 exactly nested at zero shifts', float(er@er), null['chi2'])

    # Independent central differences at the selected actual H1 endpoint.
    p = np.array(alt['parameters'])
    _, jac, _ = m1.evaluate(p)
    families = [['logN_'+str(i) for i in range(26)], ['v_'+str(i) for i in range(26)], ['logb_'+str(i) for i in range(26)]]
    probes = [max(labels, key=lambda name: np.linalg.norm(jac[:, m1.labels.index(name)])) for labels in families]
    probes += ['logb_14', 'logb_16', 'shift_1608', 'shift_2382']
    probes = list(dict.fromkeys(probes))
    derivatives = []
    for label in probes:
        j = m1.labels.index(label)
        step = .001 if label.startswith(('v_', 'shift_')) else 1e-5
        pp, pm = p.copy(), p.copy()
        pp[j] += step
        pm[j] -= step
        fd = (m1.fun(pp)-m1.fun(pm))/(2*step)
        rel = float(np.linalg.norm(fd-jac[:, j])/max(np.linalg.norm(jac[:, j]), 1e-30))
        derivatives.append({'parameter': label, 'step': step, 'relative_L2_error': rel})
        check('endpoint finite difference '+label, rel < 2e-5, f'relative L2 {rel:.8g}')
    prose = (ROOT/'reports/archive_complete_j2331_failure_pair_cn.md').read_text()
    check('Chinese report preserves unestimated precision D', '原主分析中的精密D测量仍为未估计' in prose)
    check('Chinese report excludes time and alpha promotion', '不进入宇宙年龄曲线、不转换成α' in prose)
    check('Chinese report retains local optimum limitation', '成功ftol/xtol停止不保证全局最优' in prose)
    check('no precision profile files generated', not any('profile2382' in p.name for p in fitpaths))
    check('primary JSON unchanged after validation', sha(original_path) == protocol['initial_primary_null_sha256'])
    check('primary NPZ unchanged after validation', sha(original_path.with_suffix('.npz')) == protocol['initial_primary_null_npz_sha256'])
    return protocol, summary, records, checks_by_record, derivatives, hashes


def main():
    protocol, summary, records, counts, derivatives, hashes = validate()
    failures = [x for x in CHECKS if not x[1]]
    lines = ['# Independent validation: J233156−090802 bounded failure-case pair', '',
             f"Result: **{'FAIL' if failures else 'PASS'}**, {len(CHECKS)-len(failures)}/{len(CHECKS)} explicit checks. No fitting was performed.", '',
             'The frozen diagnostic is a bounded test of whether two regional displacements repair the previously failed conventional-profile quality screen. The original precision descriptor remains unestimated. The primary gate does not constitute a population-wide no-variation test.', '',
             '## Protocol and provenance', '',
             'The runner writes a non-overwritable protocol before the first additional fit. It uses one H1 fit, one H0 cross-start, and at most one H1 cross-start when H0 improves by more than 0.001. The executed branch contains two new fits; the cross-start did not improve H0 enough to trigger the optional fit. No uncertainty profile or precision datum was generated.', '',
             f"Frozen protocol SHA-256: `{sha(OUT/'protocol.json')}`.", '',
             'Current configuration, neighboring-opacity adapter, shared engine and runner hashes match the frozen protocol. The original selected H0 JSON and NPZ still match their pre-diagnostic hashes. Additional results occupy their own directory.', '',
             '## Data and model replay', '',
             'The original H0 and both new fits use 26 shared gas components, 21 subpixels per native pixel, and exactly the same 1,176 globally unique native coadd pixels: 392 in each of the 1608, 2374 and 2382 regions. Independently recomputed window indices match all saved source indices; no mask change or pixel duplication occurred. Saved flux and expected-fluctuation errors match the full coadd exactly. Same-gas 1611 opacity remains inside the 1608 region without a separate likelihood.', '',
             'Every saved model flux, residual, full weighted Jacobian, objective, per-line objective, nominal degrees of freedom, AICc, parameter bound and source hash was replayed. Initial parameters reproduce the documented preceding fit after the same bounded projection. H0 is exactly nested in H1 at zero regional shifts. Independent direct gradients reproduce the reported first-order diagnostics.', '',
             '| Record | χ² | Nominal dof | χ²/dof | Successful stop | Explicit replay checks |',
             '|---|---:|---:|---:|---|---:|']
    for name, r in records.items():
        lines.append(f"| {name} | {r['chi2']:.9f} | {r['nominal_ndf']} | {r['chi2_per_ndf']:.9f} | {r['optimizer_success']} | {counts[name]} |")
    lines += ['', f"Selected conditional Δχ² = **{summary['delta_chi2']:.9f}**, for two extra regional shifts. The original H0 remains the better of the two H0 endpoints. H1 still fails the same global and per-line numeric quality thresholds. This limited diagnostic does not repair the overall mismatch within the frozen architecture and ±1 km/s displacement bounds.", '',
              '## Independent endpoint derivative probes', '',
              'Central differences were evaluated at the actual selected H1 endpoint for the strongest weighted-Jacobian direction in each gas-parameter family, the two active Doppler-width directions, and both regional shifts. Symmetric perturbations evaluate the smooth forward model even when one side is infinitesimally outside an optimization bound; no optimizer was called.', '',
              '| Parameter | Central step | Relative L2 error |', '|---|---:|---:|']
    for x in derivatives:
        lines.append(f"| {x['parameter']} | {x['step']:.6g} | {x['relative_L2_error']:.8g} |")
    lines += ['', '## Interpretation limitations', '',
              'Successful ftol termination does not certify a global optimum or negligible first-order gradients. Active gas-width and zero-level boundaries remain present, and the saved gradient diagnostics remain nonzero. The 1608 nuisance moves a region containing both 1608 and 1611; it is not an isolated-transition frequency measurement. The target 2382−2374 offset remains a nonprecision diagnostic. No Gaussian discovery significance, alpha value, cosmic-age datum, or logarithmic-time fit is supported by this validation.', '',
              'The finite search, chosen gas architecture, diagonal expected errors, wavelength calibration, line spread functions and blends remain limitations. The quality-screen failure itself is not evidence that every screened absorber has no physical variation.', '',
              '## Exact artifact hashes', '', '| Artifact | SHA-256 |', '|---|---|']
    for path, digest in hashes.items():
        lines.append(f'| `{path}` | `{digest}` |')
    lines += ['', '## Failed checks', '']
    lines += [f'- {name}: {detail}' for name, _, detail in failures] or ['None.']
    report = ROOT/'reports/archive_complete_j2331_failure_pair_validation.md'
    report.write_text('\n'.join(lines)+'\n')
    print(json.dumps({'checks': len(CHECKS), 'passed': len(CHECKS)-len(failures), 'failures': failures, 'derivatives': derivatives, 'report': str(report)}, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
