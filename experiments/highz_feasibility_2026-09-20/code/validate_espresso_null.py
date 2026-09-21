#!/usr/bin/env python3
"""Independent numerical checks for the scoped ESPRESSO conventional null.

Does not rerun nonlinear fits and does not turn conditional chi-square values
into discovery claims. Uses CGS cross sections and scipy's normalized
voigt_profile, independent convolution, directional finite differences, and
finer pixel quadrature to check the implementation.
"""
from pathlib import Path
import hashlib
import json
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
from scipy.special import voigt_profile
from scipy.signal import fftconvolve
import espresso_conventional_null as e

ROOT = Path(__file__).resolve().parents[1]
checks = []

def check(name, passed, **details):
    checks.append(dict(name=name, passed=bool(passed), **details))
    print(('PASS ' if passed else 'FAIL ') + name, flush=True)

def normalized_difference(a, b):
    return float(np.linalg.norm(a-b) / max(np.linalg.norm(b), 1e-15))

def independent_profile(model, p, L, ki):
    # Gaussian-CGS charge, electron mass, speed of light; Angstrom to cm.
    charge=4.803204712570263e-10
    mass=9.1093837139e-28
    light=2.99792458e10
    coefficient=np.pi*charge**2/(mass*light)
    N, velocity, doppler=model.unpack(p)
    strength_name='log_strength_'+str(L['key'])
    opacity_scale=np.exp(p[model.labels.index(strength_name)]) if strength_name in model.labels else 1.
    line_shift=(p[model.labels.index('shift_'+str(L['key']))]
                if L['key'] in model.shift_keys else 0.)
    wavelength=L['ref']*(1+e.ZREF)*np.exp(L['grid']/e.C)*1e-8
    one_plus_z=(1+e.ZREF)*np.exp((velocity+line_shift)/e.C)
    frequency=light*one_plus_z[None,:]/wavelength[:,None]
    optical_depth=np.zeros(len(wavelength))
    for lam, strength, gamma, atom_mass, q in L['atom']:
        center=light/(lam*1e-8)
        sigma=doppler*1e5/(lam*1e-8)/np.sqrt(2.)
        line_density=voigt_profile(frequency-center,sigma[None,:],gamma/(4*np.pi))
        optical_depth += (coefficient*strength*opacity_scale*N[None,:]*line_density).sum(axis=1)
    intrinsic=np.exp(-optical_depth)
    width=L['fwhm']/np.sqrt(8*np.log(2))/(model.dv/model.oversample)
    radius=int(6*width+.5)
    offsets=np.arange(-radius,radius+1)
    kernel=np.exp(-.5*(offsets/width)**2)
    kernel /= kernel.sum()
    convolved=fftconvolve(np.pad(intrinsic,radius,mode='edge'),kernel,mode='valid')
    integrated=convolved.reshape(-1,model.oversample).mean(axis=1)[L['pad']:-L['pad']]
    c0,c1=p[model.ngas+2*ki:model.ngas+2*ki+2] if model.continuum else (0.,0.)
    fixed_zero={2260:0.,2344:.002,2374:0.,2382:.002,2586:0.,2600:.003}[L['key']]
    zero_name='zero_adjustment_'+str(L['key'])
    if zero_name in model.labels:fixed_zero+=p[model.labels.index(zero_name)]
    fixed_continuum={2260:1.,2344:1.,2374:.994,2382:1.,2586:.996,2600:1.}[L['key']]
    return (fixed_zero+(1-fixed_zero)*integrated)*fixed_continuum*(1+c0+c1*L['x'])

def main():
    model=e.Model(free_shifts=True)
    p,lo,hi=model.initial()
    # Exercise derivatives away from zero continuum and zero shift.
    p[model.ngas:model.ngas+model.ncont]=np.linspace(-.001,.001,model.ncont)
    p[-len(model.shift_keys):]=np.array([.017,-.032,.008,.035,-.019])
    residual,J,profiles=model.evaluate(p)
    constant=np.sqrt(np.pi)*4.803204712570263e-10**2/(9.1093837139e-28*2.99792458e10)*1e-13
    check('CGS optical-depth coefficient',abs(e.TAU_CONSTANT/constant-1)<2e-8,
          model=e.TAU_CONSTANT,independent=constant,relative_difference=float(e.TAU_CONSTANT/constant-1))
    check('all initial arrays finite',np.isfinite(residual).all() and np.isfinite(J).all())
    check('shared gas dimension',model.ngas==135 and model.npar==152,
          gas_parameters=model.ngas,total_parameters=model.npar)
    check('anchored differential shifts',model.shift_keys==[2260,2344,2382,2586,2600])
    for ki,L in enumerate(model.lines):
        grid_centers=L['grid'].reshape(-1,model.oversample).mean(axis=1)[L['pad']:-L['pad']]
        offset=np.max(np.abs(grid_centers-L['v']))
        check(f"pixel centering {L['key']}",offset<1e-7,max_error_km_s=float(offset))
        independent=independent_profile(model,p,L,ki)
        error=np.max(np.abs(independent-profiles[ki]))
        check(f"independent CGS/Voigt/convolution {L['key']}",error<2e-7,max_flux_difference=float(error))
        check(f"finite positive noise on accepted pixels {L['key']}",
              np.isfinite(L['error'][L['good']]).all() and (L['error'][L['good']]>0).all(),
              pixels=int(L['good'].sum()))
    selected=[0,13,37,44,45,58,82,89,90,103,127,134,135,136,145,146,147,151]
    for index in selected:
        step=1e-5 if model.labels[index].startswith(('logN','continuum')) else 1e-4
        plus=p.copy();minus=p.copy();plus[index]+=step;minus[index]-=step
        finite=(model.fun(plus)-model.fun(minus))/(2*step)
        rel=normalized_difference(finite,J[:,index])
        check('Jacobian '+model.labels[index],rel<6e-5,relative_l2_error=rel,step=step)
    # Random mixed directions exercise column ordering and gas/line coupling.
    rng=np.random.default_rng(19092026)
    for trial in range(3):
        direction=rng.normal(size=model.npar)
        direction[model.ngas:model.ngas+model.ncont]*=.01
        direction/=np.linalg.norm(direction)
        step=1e-4
        numerical=(model.fun(p+step*direction)-model.fun(p-step*direction))/(2*step)
        expected=J@direction
        rel=normalized_difference(numerical,expected)
        check(f'Jacobian mixed direction {trial}',rel<3e-5,relative_l2_error=rel)
    # Use true pixel uncertainties to express quadrature discrepancies.
    fine_results={}
    for sampling in [18,27]:
        fine=e.Model(free_shifts=True,oversample=sampling)
        rf,Jf,pf=fine.evaluate(p)
        delta=rf-residual
        fine_results[sampling]=dict(residual=rf,profiles=pf)
        check(f'oversampling 9 versus {sampling}',float(delta@delta)<.01,
              noise_weighted_squared_model_difference=float(delta@delta),
              max_difference_in_sigma=float(np.max(np.abs(delta))),
              chi2_difference=float(rf@rf-residual@residual))
    delta=fine_results[27]['residual']-fine_results[18]['residual']
    check('oversampling 18 versus 27',float(delta@delta)<.001,
          noise_weighted_squared_model_difference=float(delta@delta),
          max_difference_in_sigma=float(np.max(np.abs(delta))))
    # A shared shift has exactly the summed velocity derivative before gauge fixing.
    for L,sl in zip(model.lines,model.slices):
        key=L['key']
        if key in model.shift_keys:
            a=J[sl,model.labels.index('shift_'+str(key))]
            b=J[sl,45:90].sum(axis=1)
            check(f'common-shift derivative identity {key}',normalized_difference(a,b)<1e-12,
                  relative_l2_error=normalized_difference(a,b))
    # New conventional-control branches must also pass away from their defaults.
    extended=e.Model(free_shifts=True,strength_nuisance=True,zero_nuisance=True)
    ep,_,_=extended.initial()
    for i,label in enumerate(extended.labels):
        if label.startswith('log_strength_'):ep[i]=.013
        if label.startswith('zero_adjustment_'):ep[i]=.0013
        if label.startswith('shift_'):ep[i]=.021
    er,eJ,eprofiles=extended.evaluate(ep)
    for ki,L in enumerate(extended.lines):
        independent=independent_profile(extended,ep,L,ki)
        error=float(np.max(np.abs(independent-eprofiles[ki])))
        check(f"extended independent profile {L['key']}",error<2e-7,max_flux_difference=error)
    for label in ['log_strength_2260','log_strength_2382','log_strength_2600',
                  'zero_adjustment_2260','zero_adjustment_2374','zero_adjustment_2600',
                  'shift_2382']:
        index=extended.labels.index(label)
        step=1e-5 if not label.startswith('shift_') else 1e-4
        plus=ep.copy();minus=ep.copy();plus[index]+=step;minus[index]-=step
        finite=(extended.fun(plus)-extended.fun(minus))/(2*step)
        rel=normalized_difference(finite,eJ[:,index])
        check('extended Jacobian '+label,rel<6e-5,relative_l2_error=rel,step=step)
    out=dict(status='PASS' if all(c['passed'] for c in checks) else 'FAIL',
             checks=checks,passed=sum(c['passed'] for c in checks),count=len(checks),
             root_source_sha256=hashlib.sha256(Path(e.__file__).read_bytes()).hexdigest(),
             validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             limitation='Implementation verification at perturbed published parameters; no complete nonlinear optimizer replay or statistical discovery calibration.')
    dest=ROOT/'results/espresso_null_numerical_validation.json'
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
    if out['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
