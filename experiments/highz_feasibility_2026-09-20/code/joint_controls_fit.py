#!/usr/bin/env python3
"""Actual nonlinear fits combining frozen core masks with empirical covariance.

Uses unchanged validated core and mask implementations. Covariance is rebuilt
on retained native indices, so removing pixels never joins unrelated neighbors.
Each hypothesis pair has identical masks, covariance and nuisance freedom.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import time
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import least_squares
from espresso_saturation_controls import CoreMaskModel
from espresso_conventional_null import ROOT

OUT = ROOT / 'results/joint_controls'
CONTROL = ROOT / 'results/noise_covariance/controls.json'
BASE = ROOT / 'results'
CASES = [('core01_gls', .1, 0., 'core01'), ('core03_gls', .3, 1., 'core03_padded')]


class JointModel(CoreMaskModel):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        assert self.rho == 0
        c = json.loads(CONTROL.read_text())['primary']
        acf = np.asarray(c['acf'])
        self.kernel = acf * (1 - np.arange(len(acf)) / len(acf))
        self.control_variance = c['variance']
        self.chols = []
        for line in self.lines:
            idx = line['idx'][line['good']]
            lag = np.abs(idx[:, None] - idx[None, :])
            covariance = np.zeros(lag.shape)
            inside = lag < len(self.kernel)
            covariance[inside] = self.kernel[lag[inside]]
            self.chols.append(cholesky(covariance, lower=True, check_finite=False))
        self.white_cache = None

    def evaluate(self, p, jac=True):
        if self.white_cache is not None and np.array_equal(p, self.white_cache[0]):
            return self.white_cache[1:]
        r, J, profiles = super().evaluate(p, jac)
        rw = np.empty_like(r); Jw = np.empty_like(J)
        for sl, chol in zip(self.slices, self.chols):
            rw[sl] = solve_triangular(chol, r[sl], lower=True, check_finite=False)
            Jw[sl] = solve_triangular(chol, J[sl], lower=True, check_finite=False)
        self.white_cache = (p.copy(), rw, Jw, profiles)
        return self.white_cache[1:]


def mapped_start(model, path):
    p, lo, hi = model.initial()
    old = json.loads(Path(path).read_text())
    lookup = dict(zip(old['labels'], old['parameters']))
    return np.clip([lookup.get(k, v) for k, v in zip(model.labels, p)], lo + 1e-9, hi - 1e-9)


def fit(name, free, threshold, padding, start, max_nfev=600):
    OUT.mkdir(parents=True, exist_ok=True)
    model = JointModel(free_shifts=free, threshold=threshold, padding_km_s=padding)
    _, lo, hi = model.initial()
    p = mapped_start(model, start)
    initial = p.copy(); calls = 0; began = time.time()
    def fun(x):
        nonlocal calls
        r = model.evaluate(x)[0]; calls += 1
        if calls % 25 == 0:
            print(name, calls, float(r @ r), flush=True)
            (OUT / (name + '_checkpoint.json')).write_text(json.dumps(dict(
                labels=model.labels, parameters=x.tolist(), chi2=float(r @ r), calls=calls)))
        return r
    print('START', name, model.ndata, model.npar, str(start), flush=True)
    opt = least_squares(fun, p, jac=model.jac, bounds=(lo, hi), x_scale='jac',
                        ftol=1e-7, xtol=1e-8, gtol=1e-6, max_nfev=max_nfev)
    r, J, profiles = model.evaluate(opt.x)
    singular = np.linalg.svd(J, compute_uv=False)
    dependencies = [CONTROL, BASE/'espresso_null/null_cross.npz', Path(__file__),
                    Path(__file__).with_name('espresso_saturation_controls.py'),
                    Path(__file__).with_name('espresso_conventional_null.py')]
    result = dict(name=name, free_shifts=free, threshold=threshold, padding_km_s=padding,
                  ndata=model.ndata, npar=model.npar, chi2=float(r @ r),
                  optimizer_success=bool(opt.success), status=int(opt.status), message=opt.message,
                  nfev=opt.nfev, optimality=float(opt.optimality), elapsed_seconds=time.time()-began,
                  labels=model.labels, parameters=opt.x.tolist(), initial_parameters=initial.tolist(),
                  start_file=str(Path(start).relative_to(ROOT)), mask_audit=model.mask_audit,
                  covariance_kernel=model.kernel.tolist(), continuum_variance=model.control_variance,
                  shifts_m_s={str(k):1000*float(opt.x[model.labels.index(f'shift_{k}')]) for k in model.shift_keys},
                  jacobian_rank=int(sum(singular > singular[0]*1e-10)),
                  active_bounds=[model.labels[i] for i in range(model.npar)
                                 if min(opt.x[i]-lo[i], hi[i]-opt.x[i]) < 1e-5],
                  source_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in dependencies},
                  interpretation='Exploratory paired nonlinear fit applying BOTH a frozen strong-core mask and the continuum-informed tapered covariance. Local constrained optimization, conditional covariance transfer, no calibrated discovery significance or cosmic-time-law inference.')
    (OUT / (name + '.json')).write_text(json.dumps(result, indent=2)+'\n')
    arrays = dict(parameters=opt.x, residuals=r, jacobian=J)
    for line, profile in zip(model.lines, profiles):
        for key in ['v','wave','flux','error','good','idx']:
            arrays[f'{line["key"]}_{key}'] = line[key]
        arrays[f'{line["key"]}_model'] = profile
    np.savez_compressed(OUT / (name + '.npz'), **arrays)
    print('DONE', name, result['chi2'], result['optimizer_success'], flush=True)
    return result


def run_case(case):
    name, threshold, padding, previous = case
    old = json.loads((BASE/'saturation_controls'/f'{previous}_comparison.json').read_text())
    trials = []
    def run(suffix, free, start):
        value = fit(name+'_'+suffix, free, threshold, padding, start)
        trials.append(value)
        if not value['optimizer_success']:
            value = fit(name+'_'+suffix+'_continued', free, threshold, padding, OUT/(value['name']+'.json'))
            trials.append(value)
        return value
    # Choose the lower INITIAL objective from two transparently recorded sources;
    # later opposite-hypothesis cross-starts provide a separate local basin check.
    candidates = [BASE/'saturation_controls'/(old['null']+'.json'), BASE/'empirical_gls/tapered_null_cross.json']
    model = JointModel(threshold=threshold, padding_km_s=padding)
    initial_scores = []
    for source in candidates:
        r = model.evaluate(mapped_start(model, source))[0]
        initial_scores.append(dict(file=str(source.relative_to(ROOT)), chi2=float(r @ r)))
    start = candidates[int(np.argmin([v['chi2'] for v in initial_scores]))]
    a = run('null', False, start)
    b = run('alternative', True, OUT/(a['name']+'.json'))
    c = run('alternative_cross', True, BASE/'saturation_controls'/(old['alternative']+'.json'))
    best_alt = min([r for r in trials if r['free_shifts']], key=lambda r:r['chi2'])
    run('null_cross', False, OUT/(best_alt['name']+'.json'))
    null = min([r for r in trials if not r['free_shifts']], key=lambda r:r['chi2'])
    alt = min([r for r in trials if r['free_shifts']], key=lambda r:r['chi2'])
    summary = dict(case=name, null=null['name'], alternative=alt['name'],
                   ndata=null['ndata'], added_parameters=alt['npar']-null['npar'],
                   delta_chi2=null['chi2']-alt['chi2'], null_success=null['optimizer_success'],
                   alternative_success=alt['optimizer_success'], shifts_m_s=alt['shifts_m_s'],
                   selection='Minimum recorded objective for each hypothesis, including unfinished endpoints if lower. Never silently omit lower capped endpoints.',
                   initial_null_candidates=initial_scores, mask_audit=null['mask_audit'],
                   recorded_trials=[{k:r[k] for k in ['name','free_shifts','chi2','optimizer_success','nfev','optimality']} for r in trials])
    (OUT/(name+'_comparison.json')).write_text(json.dumps(summary, indent=2)+'\n')
    return summary


if __name__ == '__main__':
    with ProcessPoolExecutor(max_workers=2) as pool:
        for comparison in pool.map(run_case, CASES):
            print(json.dumps(comparison), flush=True)
