#!/usr/bin/env python3
"""Independent implementation and saved-array audit for the UVES comparison.

Does not assert nonlinear convergence, calibrated significance, or a physical
change. Read-only on the fit products except for the adapter's metadata refresh.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import json, hashlib, tempfile
import numpy as np
from scipy.special import voigt_profile
from scipy.signal import fftconvolve
import cross_instrument_uves as u

ROOT=u.ROOT
OUT=ROOT/'results/cross_instrument'
CHECKS=[]

def check(name,passed,**kw):
    CHECKS.append(dict(name=name,passed=bool(passed),**kw))
    print(('PASS ' if passed else 'FAIL ')+name,flush=True)

def rel(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-14))

def independent_profile(model,p,line,ki):
    # Independent frequency-domain normalized Voigt and FFT convolution.
    charge=4.803204712570263e-10;electron_mass=9.1093837139e-28;light=2.99792458e10
    cross_section=np.pi*charge**2/(electron_mass*light)
    n,velocity,doppler=model.unpack(p)
    key=line['key'];label='shift_'+str(key)
    delta=p[model.labels.index(label)] if label in model.labels else 0.
    opacity_label='log_strength_'+str(key)
    opacity=np.exp(p[model.labels.index(opacity_label)]) if opacity_label in model.labels else 1.
    wavelength=line['ref']*(1+u.base.ZREF)*np.exp(line['grid']/u.base.C)*1e-8
    frequency=light*(1+u.base.ZREF)*np.exp((velocity+delta)[None,:]/u.base.C)/wavelength[:,None]
    tau=np.zeros(len(wavelength))
    for lam,f,gamma,mass,q in line['atom']:
        center=light/(lam*1e-8)
        sigma=doppler*1e5/(lam*1e-8)/np.sqrt(2.)
        tau+=(cross_section*f*opacity*n[None,:]*voigt_profile(frequency-center,sigma[None,:],gamma/(4*np.pi))).sum(axis=1)
    transmission=np.exp(-tau)
    width=line['fwhm']/np.sqrt(8*np.log(2))/(model.dv/model.oversample)
    radius=int(6*width+.5);x=np.arange(-radius,radius+1)
    kernel=np.exp(-.5*(x/width)**2);kernel/=kernel.sum()
    conv=fftconvolve(np.pad(transmission,radius,mode='edge'),kernel,mode='valid')
    pixel=conv.reshape(-1,model.oversample).mean(axis=1)[line['pad']:-line['pad']]
    c0,c1=p[model.ngas+2*ki:model.ngas+2*ki+2]
    zero=p[model.labels.index('zero_adjustment_'+str(key))]
    return (zero+(1-zero)*pixel)*(1+c0+c1*line['x'])

def construct(free=True,flex=True,oversample=7,rho=0.,**extra):
    cls=u.FlexibleLSFModel if flex else u.ORIGINAL_MODEL
    return cls(free_shifts=free,oversample=oversample,exclude=[2260,2344,2586],rho=rho,
               zero_nuisance=True,velocity_span=5.,continuum_span=.03,**extra)

def main():
    # The adapter writes metadata on load; isolate that write from active fits.
    with tempfile.TemporaryDirectory(prefix='uves_validation_metadata_') as tmp:
        u.OUT=Path(tmp)
        u.configure()
    model=construct();p,lo,hi=model.initial()
    for i,label in enumerate(model.labels):
        if label.startswith('continuum'):p[i]=.001*(1 if i%2 else -1)
        elif label.startswith('zero_adjustment'):p[i]=.0013
        elif label.startswith('shift_'):p[i]=.039
        elif label.startswith('log_lsf'):p[i]=-.15+.02*(i-model.base_npar)
    r,j,profiles=model.evaluate(p)
    check('finite initialized arrays',np.isfinite(r).all() and np.isfinite(j).all())
    check('three selected lines and anchor',model.keys==[2374,2382,2600] and model.shift_keys==[2382,2600])
    check('data and nuisance dimensions',model.ndata==511 and model.npar==149 and model.ngas==135,
          ndata=model.ndata,npar=model.npar)
    null=construct(free=False)
    check('nested parameter parity',set(null.labels)==set(model.labels)-{'shift_2382','shift_2600'})
    p_null,_,_=null.initial();lookup=dict(zip(model.labels,p))
    p_null=np.array([lookup[k] for k in null.labels]);p_zero=p.copy()
    for k in model.shift_keys:p_zero[model.labels.index(f'shift_{k}')]=0.
    check('nested residual parity at zero shifts',np.array_equal(null.evaluate(p_null)[0],model.evaluate(p_zero)[0]))
    # Restore widths, profiles, and Jacobian after the nesting call.
    r,j,profiles=model.evaluate(p)
    src=np.load(ROOT/'data/processed/J051707-441055_squad_dr1.npz')
    for ki,line in enumerate(model.lines):
        key=line['key'];idx=line['idx'];good=line['good']
        check(f'source flux and error {key}',np.array_equal(line['flux'],src['flux_normalized'][idx]) and
              np.array_equal(line['error'],src['error_normalized'][idx]))
        expected_good=src['valid'][idx]&np.isfinite(src['flux_normalized'][idx])&(src['error_normalized'][idx]>0)
        check(f'unchanged clean-line mask {key}',np.array_equal(good,expected_good),valid_pixels=int(good.sum()))
        centers=line['grid'].reshape(-1,model.oversample).mean(axis=1)[line['pad']:-line['pad']]
        error=float(np.max(np.abs(centers-line['v'])))
        check(f'pixel quadrature centering {key}',error<1e-7,max_error_km_s=error)
        error=float(np.max(np.abs(profiles[ki]-independent_profile(model,p,line,ki))))
        check(f'independent frequency Voigt and FFT profile {key}',error<2e-7,max_flux_error=error)
    indices=[0,13,37,44,45,58,82,89,90,103,127,134]+list(range(135,model.npar))
    for idx in indices:
        label=model.labels[idx];step=1e-5 if label.startswith(('logN','logb','continuum','zero_')) else 1e-4
        plus=p.copy();minus=p.copy();plus[idx]+=step;minus[idx]-=step
        fd=(model.evaluate(plus)[0]-model.evaluate(minus)[0])/(2*step)
        error=rel(fd,j[:,idx]);check('derivative '+label,error<7e-5,relative_l2_error=error,step=step)
    rng=np.random.default_rng(20260921)
    for trial in range(4):
        d=rng.normal(size=model.npar);d[135:144]*=.01;d/=np.linalg.norm(d);step=1e-4
        fd=(model.evaluate(p+step*d)[0]-model.evaluate(p-step*d)[0])/(2*step)
        error=rel(fd,j@d);check(f'mixed derivative {trial}',error<7e-5,relative_l2_error=error)
    # AR(1) branches must whiten both residuals and analytic/width derivatives.
    corr=construct(rho=.3);cr,cj,_=corr.evaluate(p)
    d=rng.normal(size=model.npar);d[135:144]*=.01;d/=np.linalg.norm(d);step=1e-4
    fd=(corr.evaluate(p+step*d)[0]-corr.evaluate(p-step*d)[0])/(2*step)
    error=rel(fd,cj@d);check('AR(1) mixed derivative',error<7e-5,relative_l2_error=error)
    extended=construct(strength_nuisance=True,error_row=2)
    ep,_,_=extended.initial()
    for i,label in enumerate(extended.labels):
        if label.startswith('log_strength'):ep[i]=.021
        elif label.startswith('log_lsf'):ep[i]=-.15
        elif label.startswith('zero_adjustment'):ep[i]=.001
        elif label.startswith('shift_'):ep[i]=.027
    er,ej,em=extended.evaluate(ep)
    for ki,line in enumerate(extended.lines):
        key=line['key'];error=float(np.max(np.abs(em[ki]-independent_profile(extended,ep,line,ki))))
        check(f'opacity independent profile {key}',error<2e-7,max_flux_error=error)
        check(f'expected-fluctuation errors {key}',np.array_equal(line['error'],src['expected_fluctuation_normalized'][line['idx']]))
    for label in ['log_strength_2382','log_strength_2600','log_lsf_scale_2382','shift_2382']:
        i=extended.labels.index(label);step=1e-4;plus=ep.copy();minus=ep.copy();plus[i]+=step;minus[i]-=step
        fd=(extended.evaluate(plus)[0]-extended.evaluate(minus)[0])/(2*step)
        error=rel(fd,ej[:,i]);check('opacity derivative '+label,error<7e-5,relative_l2_error=error)
    # Recompute each finalized result; ignore mutable checkpoints.
    audited=[]
    for path in sorted(OUT.glob('uves_*.json')):
        if path.name.endswith('_checkpoint.json'):continue
        data=json.loads(path.read_text())
        if 'parameters' not in data or 'configuration' not in data:continue
        if not path.with_suffix('.npz').exists():continue
        flex=any(x.startswith('log_lsf_scale') for x in data['labels'])
        cls=u.FlexibleLSFModel if flex else u.ORIGINAL_MODEL
        saved=cls(free_shifts=data['free_shifts'],**data['configuration'])
        sp=np.array(data['parameters']);sr,sj,sm=saved.evaluate(sp)
        _,lower,upper=saved.initial()
        check(path.stem+' labels',saved.labels==data['labels'])
        check(path.stem+' parameter bounds',np.all(sp>=lower-1e-10) and np.all(sp<=upper+1e-10))
        check(path.stem+' accounting',saved.ndata==data['ndata'] and saved.npar==data['npar'] and
              data['nominal_ndf']==data['ndata']-data['npar'])
        check(path.stem+' chi2',abs(float(sr@sr)-data['chi2'])<1e-7,recomputed=float(sr@sr),saved=data['chi2'])
        z=np.load(path.with_suffix('.npz'))
        check(path.stem+' residual arrays',np.max(np.abs(sr-z['residuals']))<1e-8)
        check(path.stem+' Jacobian arrays',np.max(np.abs(sj-z['jacobian']))<1e-7)
        for ki,(line,prof) in enumerate(zip(saved.lines,sm)):
            key=line['key'];check(path.stem+f' profile {key}',np.max(np.abs(prof-z[f'{key}_model']))<1e-10)
            check(path.stem+f' saved data arrays {key}',all(np.array_equal(line[k],z[f'{key}_{k}']) for k in ['v','wave','flux','error','good']))
            if data['configuration']['oversample']>=21:
                independent=independent_profile(saved,sp,line,ki)
                error=float(np.max(np.abs(prof-independent)))
                check(path.stem+f' independent fitted CGS profile {key}',error<2e-7,max_flux_error=error)
        # Numerical quadrature assessment at actual fitted narrow components.
        if flex:
            sampling=data['configuration']['oversample']
            kw=dict(data['configuration']);kw['oversample']=max(21,sampling*2+7)
            fine=cls(free_shifts=data['free_shifts'],**kw);fr,_,_=fine.evaluate(sp)
            delta=fr-sr
            # Products at 7 samples are retained as exploratory records; their
            # refinement difference is measured, never represented as convergence.
            label=f" pixel quadrature {sampling}-to-{kw['oversample']}"
            accepted=float(delta@delta)<.05 if sampling>=21 else np.isfinite(delta).all()
            check(path.stem+label,accepted,
                  squared_noise_weighted_difference=float(delta@delta),
                  max_sigma_difference=float(np.max(np.abs(delta))),chi2_difference=float(fr@fr-sr@sr),
                  assessment='refinement convergence check' if sampling>=21 else 'coarse exploratory characterization only; use >=21 for final comparison')
        audited.append(dict(name=path.stem,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            optimizer_success=data['optimizer_success'],optimality=data['optimality'],
                            chi2=data['chi2'],npar=data['npar'],rank=data['jacobian_rank'],
                            oversample=data['configuration']['oversample'],
                            active_bounds=data['active_bounds'],shifts=data['shifts_km_s']))
    out=dict(status='PASS' if all(x['passed'] for x in CHECKS) else 'FAIL',
             passed=sum(x['passed'] for x in CHECKS),count=len(CHECKS),checks=CHECKS,
             audited_results=audited,adapter_sha256=hashlib.sha256(Path(u.__file__).read_bytes()).hexdigest(),
             validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             inference_status='Not validated: implementation checks do not establish optimizer convergence or calibrated statistical inference.',
             interpretation='Implementation/array verification only. Optimizer failures, flexible/bound-active LSF, incomplete calibration-provenance reproduction and weakly identified or bound-active gas parameters prevent calibrated physical inference.')
    (OUT/'validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['checks','audited_results']},indent=2))
    if out['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
