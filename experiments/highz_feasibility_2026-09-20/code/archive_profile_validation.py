#!/usr/bin/env python3
"""Independent implementation audit of the exploratory archive profile pilot.

Does not fit data or certify the physical adequacy of the chosen gas model.
Reconstructs native data from raw FITS, finite-differences all parameter types,
uses scipy's separately implemented Voigt profile and an explicit Gaussian
convolution for the forward profile, and replays stored fits and quadrature.
"""
import hashlib
import io
import json
import math
import sys
import tarfile
from pathlib import Path

import numpy as np
from scipy.special import voigt_profile

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / '.deps'))
from astropy.io import fits
import archive_profile_pilot as pilot

OUT = BASE / 'results/archive_expansion/profile_pilot'
CHECKS = []
QUADRATURE = []
DERIVATIVES = []


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check(name, passed, detail=None):
    CHECKS.append(dict(name=name, passed=bool(passed), detail=detail))


def close(name, actual, expected, rtol=2e-11, atol=2e-11):
    a, b = np.asarray(actual), np.asarray(expected)
    ok = a.shape == b.shape and np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True)
    detail = float(np.nanmax(np.abs(a.astype(float) - b.astype(float)))) if a.shape == b.shape and a.size else None
    check(name, ok, detail)


def independent_atoms(key):
    rows = []
    for line in pilot.ATOMIC.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 6 and fields[0] == 'FeII' and key <= float(fields[1]) < key + 1:
            rows.append([float(x) for x in fields[1:6]])
    return np.asarray(rows)


def independent_profile(model, p, line):
    # voigt_profile uses a normalized Gaussian sigma, whereas the source's
    # Re[wofz] convention has sigma=1/sqrt(2); multiplying sqrt(pi) converts it.
    n = model.ncomp
    velocity_grid = line['grid']
    tau = np.zeros(len(velocity_grid))
    shift = p[model.labels.index(f'shift_{line["key"]}')] if line['key'] in model.shiftkeys else 0.
    for column, velocity, log_b in zip(p[:n], p[n:2*n], p[2*n:3*n]):
        b = math.exp(log_b)
        for wave, oscillator, damping, _mass, _q in independent_atoms(line['key']):
            centre = velocity + pilot.C * math.log(wave / line['ref']) + shift
            u = -pilot.C * np.expm1((centre - velocity_grid) / pilot.C) / b
            a = damping * wave * 1e-13 / (4 * math.pi * b)
            tau += pilot.TAU * 10.**column * oscillator * wave / b * math.sqrt(math.pi) * voigt_profile(u, 1 / math.sqrt(2), a)
    offset = 3*n + 4*model.keys.index(line['key'])
    c0, c1, zero, log_fwhm = p[offset:offset+4]
    sigma = math.exp(log_fwhm) / 2.354820045 / (line['dv'] / model.oversample)
    radius = int(6. * sigma + .5)
    xx = np.arange(-radius, radius + 1)
    kernel = np.exp(-.5 * (xx / sigma)**2)
    kernel /= kernel.sum()
    fine = np.convolve(np.pad(np.exp(-tau), (radius, radius), mode='edge'), kernel, mode='valid')
    coarse = fine.reshape(-1, model.oversample).mean(axis=1)[line['pad']:-line['pad']]
    return (1 + c0 + c1*line['x']) * (zero + (1-zero)*coarse)


def audit_derivatives(model, p, label):
    residual, jacobian, profiles = model.evaluate(p)
    original = p.copy()
    p[0] += .01
    check(label + ' cache parameter copy', not np.array_equal(model.cache[0], p))
    p[:] = original
    check(label + ' cache hit same arrays', model.evaluate(p)[0] is residual)
    for j, name in enumerate(model.labels):
        step = 1e-3 if name.startswith(('v_', 'shift_')) else (1e-5 if name.startswith(('cont', 'zero')) else 1e-4)
        plus, minus = p.copy(), p.copy()
        plus[j] += step
        minus[j] -= step
        difference = (model.fun(plus) - model.fun(minus)) / (2*step)
        error = np.linalg.norm(difference-jacobian[:,j]) / max(np.linalg.norm(jacobian[:,j]), 1e-4)
        DERIVATIVES.append(dict(case=label, parameter=name, step=step, relative_l2_error=float(error)))
        check(label + ' derivative ' + name, error < 8e-5, float(error))
    for line, expected in zip(model.lines, profiles):
        actual = independent_profile(model, p, line)
        close(label + ' independent Voigt Gaussian profile ' + str(line['key']), actual, expected, atol=1e-9, rtol=1e-9)


def raw_arrays(target):
    meta = json.loads((pilot.DATA / f'{target}_metadata.json').read_text())
    path = BASE / meta['archive_file']
    check(target + ' raw tar SHA256', digest(path) == meta['archive_sha256'])
    with tarfile.open(path) as archive:
        data = archive.extractfile(f'{target}/{target}.fits').read()
    check(target + ' raw FITS SHA256', hashlib.sha256(data).hexdigest() == meta['fits_sha256'])
    with fits.open(io.BytesIO(data)) as hdus:
        arr = hdus[0].data.copy()
        header = hdus[0].header.copy()
    wave = 10.**(header['CRVAL1'] + (np.arange(arr.shape[1])+1-header['CRPIX1'])*header['CD1_1'])
    valid = (arr[4] == 1) & np.isfinite(arr[0]) & np.isfinite(arr[1]) & (arr[1] > 0)
    return wave, arr, valid


def main():
    # Direct reconstruction of the full source data and exact frozen windows.
    expected_f = {1608:.0577, 2260:.00244, 2344:.114, 2374:.0313, 2382:.32, 2586:.0691}
    for target in pilot.CONFIG:
        wave, raw, valid = raw_arrays(target)
        model = pilot.Model(target, 3, True)
        for line in model.lines:
            key = line['key']; prefix = f'{target} FeII{key}'
            atom = independent_atoms(key)
            close(prefix + ' isotope rows', line['atom'], atom)
            close(prefix + ' total oscillator strength', atom[:,1].sum(), expected_f[key], atol=2e-8)
            if len(atom) == 4:
                close(prefix + ' oscillator strengths already abundance weighted', atom[:,1]/atom[:,1].sum(), [.00282,.02119,.91754,.05845], atol=2e-7)
            ref = np.average(atom[:,0], weights=atom[:,1])
            close(prefix + ' reference wavelength', line['ref'], ref)
            velocity = pilot.C*np.log(wave / (ref*(1+model.cfg['z'])))
            ids = np.flatnonzero((velocity >= model.cfg['window'][0]) & (velocity <= model.cfg['window'][1]))
            close(prefix + ' exact raw source indices', line['source_indices'], ids, rtol=0, atol=0)
            close(prefix + ' native wavelengths', line['wave'], wave[ids], atol=3e-11)
            close(prefix + ' native velocity', line['v'], velocity[ids], atol=1e-8)
            close(prefix + ' native flux', line['flux'], raw[0,ids], rtol=0, atol=0)
            close(prefix + ' native errors', line['error'], raw[1,ids], rtol=0, atol=0)
            close(prefix + ' source validity without residual clipping', line['good'], valid[ids], rtol=0, atol=0)
            check(prefix + ' consecutive native pixels', np.all(np.diff(ids)==1))
            check(prefix + ' finite positive retained errors', np.all(np.isfinite(line['error'][line['good']]) & (line['error'][line['good']]>0)))
        p, lower, upper = model.initial(0)
        # Exercise nonzero continuum, zero levels, line shifts and all widths.
        rng = np.random.default_rng(260921)
        p += rng.uniform(-.001, .001, len(p))
        for j, name in enumerate(model.labels):
            if name.startswith('shift_'): p[j] = .18
        audit_derivatives(model, p, target + ' interior n3')
        null = pilot.Model(target,3,False)
        nullp = p[:null.npar]
        zerop = p.copy(); zerop[null.npar:] = 0
        close(target + ' H0 nested exactly inside H1', null.fun(nullp), model.fun(zerop), rtol=0, atol=0)
        for a,b in zip(null.lines,model.lines):
            close(target + ' identical H0 H1 retained pixels ' + str(a['key']),a['source_indices'][a['good']],b['source_indices'][b['good']],rtol=0,atol=0)

    fit_files = sorted(x for x in OUT.glob('*.json') if '_n' in x.stem and ('_null_' in x.stem or '_shift_' in x.stem))
    for path in fit_files:
        obj = json.loads(path.read_text()); model = pilot.Model(obj['target'],obj['ncomp'],obj['free_shifts'],obj['oversample'])
        p = np.asarray(obj['parameters']); r,j,profiles = model.evaluate(p)
        prefix = obj['name']; saved = np.load(path.with_suffix('.npz'))
        check(prefix + ' parameter labels',obj['labels']==model.labels)
        check(prefix + ' model SHA256',obj['model_script_sha256']==digest(pilot.__file__))
        check(prefix + ' atomic SHA256',obj['atomic_sha256']==digest(pilot.ATOMIC))
        check(prefix + ' input window SHA256',all(digest(k)==v for k,v in obj['source_files'].items()))
        check(prefix + ' dimensions',model.ndata==obj['ndata'] and model.npar==obj['npar'] and model.ndata-model.npar==obj['nominal_ndf'])
        close(prefix + ' stored parameter vector',saved['parameters'],p,rtol=0,atol=0)
        close(prefix + ' residual replay',saved['residual'],r)
        close(prefix + ' Jacobian replay',saved['jacobian'],j)
        close(prefix + ' chi2 replay',obj['chi2'],r@r)
        close(prefix + ' chi2/ndf',obj['chi2_per_ndf'],(r@r)/(model.ndata-model.npar))
        expected_aicc = r@r+2*model.npar+2*model.npar*(model.npar+1)/(model.ndata-model.npar-1)
        close(prefix + ' AICc arithmetic',obj['aicc'],expected_aicc)
        lo,hi = np.asarray(obj['bounds_lower']),np.asarray(obj['bounds_upper'])
        check(prefix + ' fitted parameter feasibility',np.all(p>=lo) and np.all(p<=hi))
        for line,profile,sl,stated in zip(model.lines,profiles,model.slices,obj['per_line']):
            q=prefix+' '+str(line['key'])
            close(q+' saved flux profile',saved[f'{line["key"]}_model'],profile)
            close(q+' per line chi2',stated['chi2'],r[sl]@r[sl])
            for k in ['v','wave','flux','error','good','source_indices']:
                close(q+' saved '+k,saved[f'{line["key"]}_{k}'],line[k],rtol=0,atol=0)
        for k,shift in obj['shifts_m_s'].items():
            close(prefix+' shift units '+k,shift,1000*p[model.labels.index('shift_'+k)])

    for path in sorted(OUT.glob('*_comparison.json')):
        comparisons=json.loads(path.read_text())['comparisons']
        for obj in comparisons:
            left=json.loads((OUT/(obj['null']+'.json')).read_text())
            right=json.loads((OUT/(obj['alternative']+'.json')).read_text())
            prefix=path.stem+f' n{obj["ncomp"]}'
            close(prefix+' delta chi2 arithmetic',obj['delta_chi2'],left['chi2']-right['chi2'])
            check(prefix+' fit status fidelity',obj['both_success']==(left['optimizer_success'] and right['optimizer_success']))
            check(prefix+' nested parameter count',right['npar']-left['npar']==len(right['shifts_m_s']))
            objectives=[]
            for fit in [left,right]:
                p=np.asarray(fit['parameters']); models=[pilot.Model(fit['target'],fit['ncomp'],fit['free_shifts'],s) for s in [11,33,55]]
                residuals=[m.fun(p) for m in models]
                chi=[float(r@r) for r in residuals]
                objectives.append(chi)
                QUADRATURE.append(dict(name=fit['name'],oversamples=[11,33,55],chi2=chi,chi2_33_minus_11=chi[1]-chi[0],chi2_55_minus_33=chi[2]-chi[1],max_normalized_profile_difference_11_55=float(np.max(np.abs(residuals[0]-residuals[2]))),l2_normalized_profile_difference_11_55=float(np.linalg.norm(residuals[0]-residuals[2]))))
                check(fit['name']+' 33 to55 quadrature small objective change',abs(chi[2]-chi[1])<.1,float(chi[2]-chi[1]))
                if fit['ncomp']==6:
                    audit_derivatives(models[0],p.copy(),fit['name']+' fitted n6')
            d=np.array(objectives[0])-np.array(objectives[1])
            check(prefix+' delta chi2 11 to55 quadrature stable within 0.1',abs(d[2]-d[0])<.1,d.tolist())

    selected_metrics=[]
    for path in sorted(OUT.glob('*_selected_refined.json')):
        selection=json.loads(path.read_text())
        left=json.loads((OUT/(selection['null']+'.json')).read_text())
        right=json.loads((OUT/(selection['alternative']+'.json')).read_text())
        rows=json.loads((OUT/(selection['target']+'_comparison.json')).read_text())['comparisons']
        check(path.stem+' AICc selection arithmetic',selection['ncomp']==min(rows,key=lambda x:x['aicc_null'])['ncomp'])
        close(path.stem+' selected delta chi2',selection['delta_chi2'],left['chi2']-right['chi2'])
        check(path.stem+' selected fit status',selection['both_success']==(left['optimizer_success'] and right['optimizer_success']))
        check(path.stem+' both selected fits terminated successfully',left['optimizer_success'] and right['optimizer_success'])
        check(path.stem+' extra parameter count',selection['extra_shifts']==right['npar']-left['npar']==len(right['shifts_m_s']))
        pair=[]
        for fit in [left,right]:
            models=[pilot.Model(fit['target'],fit['ncomp'],fit['free_shifts'],s) for s in [33,55]]
            residuals=[m.fun(np.array(fit['parameters'])) for m in models]
            chi=[float(r@r) for r in residuals]
            check(fit['name']+' selected quadrature33',fit['oversample']==33)
            check(fit['name']+' selected 33 to55 objective within0.01',abs(chi[1]-chi[0])<.01,float(chi[1]-chi[0]))
            for line in models[0].lines:
                close(fit['name']+' selected retained native pixels '+str(line['key']),line['source_indices'][line['good']],models[1].lines[models[0].keys.index(line['key'])]['source_indices'][line['good']],rtol=0,atol=0)
            pair.append(dict(name=fit['name'],chi2_33=chi[0],chi2_55=chi[1],delta_chi2_55_minus33=chi[1]-chi[0],maximum_profile_difference_in_error_units=float(np.max(np.abs(residuals[1]-residuals[0])))))
        selected_metrics.append(dict(target=selection['target'],ncomp=selection['ncomp'],null=pair[0],alternative=pair[1],delta_chi2_33=pair[0]['chi2_33']-pair[1]['chi2_33'],delta_chi2_55=pair[0]['chi2_55']-pair[1]['chi2_55']))
    check('both authoritative selected refined targets present',len(selected_metrics)==2)

    result=dict(scope='Independent raw-FITS, derivative, forward-profile, saved-fit and quadrature implementation audit; not physical model adequacy or global optimization certification.',passed=all(c['passed'] for c in CHECKS),n_checks=len(CHECKS),n_saved_fits=len(fit_files),model_sha256=digest(pilot.__file__),validator_sha256=digest(__file__),max_derivative_relative_l2_error=max(x['relative_l2_error'] for x in DERIVATIVES),quadrature=QUADRATURE,selected_refined=selected_metrics,derivatives=DERIVATIVES,failures=[c for c in CHECKS if not c['passed']],checks=CHECKS)
    (OUT/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# Independent archive profile implementation validation','',f'Passed: **{result["passed"]}**. {len(CHECKS)} checks; {len(fit_files)} saved fits replayed.','',f'Maximum finite-difference Jacobian relative L2 error: {result["max_derivative_relative_l2_error"]:.6g}.','', 'Validation reconstructs native wavelengths, exact indices, errors and validity from the original FITS tar members. It checks isotope oscillator strengths without multiplying abundances twice, all parameter derivatives, an independently evaluated normalized Voigt profile with explicit Gaussian convolution, cache parameter copies, identical H0/H1 pixels, saved objective values, parameter feasibility, and AICc arithmetic.','', 'Quadrature comparisons evaluate the same saved parameters; they do not re-optimize at higher sampling. The complete metrics are in `results/archive_expansion/profile_pilot/independent_validation.json`.','', '| Fit | χ² at 11 | 33−11 | 55−33 | Largest profile difference / pixel error, 11 vs55 |','|---|---:|---:|---:|---:|']
    for q in QUADRATURE:
        lines.append(f'| {q["name"]} | {q["chi2"][0]:.6f} | {q["chi2_33_minus_11"]:.6f} | {q["chi2_55_minus_33"]:.6f} | {q["max_normalized_profile_difference_11_55"]:.6g} |')
    lines += ['', '## Authoritative selected refined fits', '', '| Target | Components | Δχ² at33 | Δχ² at55, same parameters |', '|---|---:|---:|---:|']
    for s in selected_metrics:
        lines.append(f'| {s["target"]} | {s["ncomp"]} | {s["delta_chi2_33"]:.9f} | {s["delta_chi2_55"]:.9f} |')
    lines += ['', 'All selected fits passed their recorded termination status checks. Higher quadrature evaluations are not additional nonlinear fits; source results at33 are the authoritative numerical optimization products.', '', 'The pilot remains exploratory: velocity windows and component counts were developed after viewing profiles; diagonal errors are used; several nuisance parameters reach bounds; instrumental and blend uncertainties remain uncalibrated. Successful scipy stopping is not proof of a global optimum. AICc is descriptive under these assumptions, not decisive physical model selection. A parameter called a relative shift is not an alpha measurement or a cosmic time trend.']
    if result['failures']:
        lines += ['', 'Failed checks:', '']+[f'- {c["name"]}: {c["detail"]}' for c in result['failures']]
    (BASE/'reports/archive_profile_validation.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['checks','derivatives','quadrature']},indent=2),flush=True)


if __name__=='__main__':
    main()
