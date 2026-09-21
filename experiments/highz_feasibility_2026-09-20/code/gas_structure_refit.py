#!/usr/bin/env python3
"""A fixed, initialization-only component-merging sensitivity architecture.

The rule is not selected using fitted line offsets. Sub-resolution components
can contain measurable shape information at high S/N; this deliberately
parsimonious architecture is a dependence test, not an assertion of truth.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import argparse, hashlib, json, time
import numpy as np
from scipy.optimize import least_squares
from espresso_conventional_null import load_sources, ROOT
from espresso_empirical_gls import GLSModel, CONTROL

OUT=ROOT/'results/gas_structure'
BASE=ROOT/'results/empirical_gls'


def merged_components(comp, threshold):
    order=np.argsort(comp[:,1], kind='stable')
    cuts=np.flatnonzero(np.diff(comp[order,1])>=threshold)+1
    groups=[g.tolist() for g in np.split(order,cuts)]
    new=[]
    for g in groups:
        sub=comp[g]; columns=10**sub[:,0]; weights=columns/columns.sum()
        velocity=float(weights@sub[:,1])
        b2=float(weights@(sub[:,2]**2+2*(sub[:,1]-velocity)**2))
        new.append([np.log10(columns.sum()),velocity,np.sqrt(b2)])
    return np.asarray(new),groups


class MergedGLSModel(GLSModel):
    def __init__(self, threshold=1.025, **kwargs):
        super().__init__(**kwargs)
        self.published_comp=self.comp.copy()
        self.threshold=threshold
        self.comp,self.groups=merged_components(self.comp,threshold)
        tail=self.labels[self.ngas:]
        self.ncomp=len(self.comp)
        self.labels=[f'logN_{i}' for i in range(self.ncomp)]
        if self.fit_velocity:self.labels += [f'v_{i}' for i in range(self.ncomp)]
        if self.fit_b:self.labels += [f'logb_{i}' for i in range(self.ncomp)]
        self.ngas=len(self.labels)
        self.labels += tail
        self.strength_offset=self.ngas+self.ncont
        self.zero_offset=self.strength_offset+len(self.strength_keys)
        self.shift_offset=self.zero_offset+len(self.zero_keys)
        self.npar=len(self.labels)
        self.cache=None;self.white_cache=None


def convert_start(model, start):
    initial,lo,hi=model.initial()
    if start is None:return initial,lo,hi
    old=json.loads(Path(start).read_text())
    lookup=dict(zip(old['labels'],old['parameters']))
    if old.get('ncomp')==model.ncomp and old.get('merge_groups')==model.groups:
        initial=np.array([lookup.get(k,val) for k,val in zip(model.labels,initial)])
    else:
        nold=sum(k.startswith('logN_') for k in old['labels'])
        if nold!=45:raise ValueError('Can only map published 45-component or identical architecture starts')
        fitted=np.array([[lookup[f'logN_{i}'],lookup[f'v_{i}'],np.exp(lookup[f'logb_{i}'])] for i in range(nold)])
        # Keep the fixed published-membership groups; never regroup by fitted centroids.
        merged=[]
        for g in model.groups:
            sub=fitted[g];N=10**sub[:,0];w=N/N.sum();v=w@sub[:,1]
            merged.append([np.log10(N.sum()),v,np.sqrt(w@(sub[:,2]**2+2*(sub[:,1]-v)**2))])
        merged=np.asarray(merged)
        for i in range(model.ncomp):
            lookup[f'logN_{i}']=merged[i,0];lookup[f'v_{i}']=merged[i,1];lookup[f'logb_{i}']=np.log(merged[i,2])
        initial=np.array([lookup.get(k,val) for k,val in zip(model.labels,initial)])
    return np.clip(initial,lo+1e-9,hi-1e-9),lo,hi


def fit(name, free=False, start=None, max_nfev=400, threshold=1.025):
    OUT.mkdir(parents=True,exist_ok=True)
    m=MergedGLSModel(free_shifts=free, threshold=threshold)
    p,lo,hi=convert_start(m,start);initial=p.copy();calls=0;began=time.time()
    def fun(p):
        nonlocal calls
        r=m.evaluate(p)[0];calls+=1
        if calls%25==0:
            print(name,calls,float(r@r),flush=True)
            checkpoint=dict(name=name,ncomp=m.ncomp,merge_groups=m.groups,labels=m.labels,parameters=p.tolist(),chi2=float(r@r))
            (OUT/(name+'_checkpoint.json')).write_text(json.dumps(checkpoint,indent=2)+'\n')
        return r
    print('START',name,m.ncomp,m.ndata,m.npar,flush=True)
    opt=least_squares(fun,p,jac=m.jac,bounds=(lo,hi),x_scale='jac',ftol=1e-7,xtol=1e-8,gtol=1e-6,max_nfev=max_nfev)
    r,J,profiles=m.evaluate(opt.x);s=np.linalg.svd(J,compute_uv=False)
    result=dict(name=name,free_shifts=free,threshold_km_s=threshold,ncomp=m.ncomp,merge_groups=m.groups,published_components=m.published_comp.tolist(),merged_initial_components=m.comp.tolist(),ndata=m.ndata,npar=m.npar,nominal_ndf=m.ndata-m.npar,chi2=float(r@r),optimizer_success=bool(opt.success),message=opt.message,nfev=opt.nfev,optimality=float(opt.optimality),elapsed_seconds=time.time()-began,labels=m.labels,parameters=opt.x.tolist(),initial_parameters=initial.tolist(),start_file=str(Path(start).relative_to(ROOT)) if start else None,jacobian_rank=int(sum(s>s[0]*1e-10)),active_bounds=[m.labels[i] for i in range(m.npar) if min(opt.x[i]-lo[i],hi[i]-opt.x[i])<1e-5],shifts_m_s={k:1000*float(opt.x[m.labels.index(f'shift_{k}')]) for k in m.shift_keys},covariance_kernel=m.kernel.tolist(),continuum_variance=m.variance,covariance_source_sha256=hashlib.sha256(CONTROL.read_bytes()).hexdigest(),core_source_sha256=hashlib.sha256((ROOT/'code/espresso_conventional_null.py').read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),interpretation='A fixed published-initialization component-merging sensitivity architecture under unchanged empirical GLS errors. Not an exhaustive gas-architecture search or a physical-discovery probability.')
    (OUT/f'{name}.json').write_text(json.dumps(result,indent=2)+'\n')
    arrays=dict(parameters=opt.x,residuals=r,jacobian=J)
    for L,profile in zip(m.lines,profiles):
        for k in ['v','wave','flux','error','good']:arrays[f'{L["key"]}_{k}']=L[k]
        arrays[f'{L["key"]}_model']=profile
    np.savez_compressed(OUT/f'{name}.npz',**arrays)
    print('DONE',name,float(r@r),bool(opt.success),opt.message,flush=True)
    return result


def continued(name, free, start=None, max_nfev=400):
    result=fit(name,free,start,max_nfev)
    if not result['optimizer_success']:
        result=fit(name+'_continued',free,OUT/(result['name']+'.json'),max_nfev)
    return result


def workflow():
    OUT.mkdir(parents=True,exist_ok=True)
    definition=MergedGLSModel()
    design=dict(rule='Connected groups in ascending published centroid order with adjacent separations < half the minimum adopted Gaussian instrumental FWHM; threshold fixed before new offsets.',threshold_km_s=definition.threshold,ncomp=definition.ncomp,groups=definition.groups,initialization='Conserve total N, N-weighted centroid and second moment b^2/2 of Doppler core. Voigt Lorentzian moments are not finite; this only specifies a Doppler-core initializer.',model_scope='Refit all merged component columns, centroids, Doppler widths, and the same 12 continuum coefficients; unchanged six transitions, pixels, atomic data, fixed LSF and empirical GLS kernel.',bounds='Same published-initializer-centered ranges as baseline: logN +/-2 dex (clipped 7..17), velocity +/-2 km/s, b 0.15..30 km/s, continuum +/-0.03. Merge membership remains fixed.',interpretation='A parsimonious dependence test; half an instrumental FWHM is not a proof that component structure is unmeasurable at high S/N. Neither independent architecture selection nor calibrated significance.',code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'design.json').write_text(json.dumps(design,indent=2)+'\n')
    a=continued('merged40_null_published',False)
    b=continued('merged40_alternative_from_null',True,OUT/(a['name']+'.json'))
    base=json.loads((BASE/'comparison.json').read_text())
    c=continued('merged40_alternative_cross',True,BASE/(base['alternative']+'.json'))
    bs=[b,c];best=min([r for r in bs if r['optimizer_success']] or bs,key=lambda r:r['chi2'])
    d=continued('merged40_null_cross',False,OUT/(best['name']+'.json'))
    aa=[a,d];an=min([r for r in aa if r['optimizer_success']] or aa,key=lambda r:r['chi2'])
    oldnull=json.loads((BASE/(base['null']+'.json')).read_text());oldalt=json.loads((BASE/(base['alternative']+'.json')).read_text())
    summary=dict(null=an['name'],alternative=best['name'],null_success=an['optimizer_success'],alternative_success=best['optimizer_success'],ncomp=an['ncomp'],ndata=an['ndata'],null_parameters=an['npar'],alternative_parameters=best['npar'],null_chi2=an['chi2'],alternative_chi2=best['chi2'],delta_chi2=an['chi2']-best['chi2'],added_parameters=best['npar']-an['npar'],shifts_m_s=best['shifts_m_s'],baseline45=dict(null=base['null'],alternative=base['alternative'],null_chi2=oldnull['chi2'],alternative_chi2=oldalt['chi2'],delta_chi2=base['delta_chi2']),all_final_attempts=[{k:r[k] for k in ['name','chi2','optimizer_success','nfev','optimality']} for r in [a,b,c,d]],limitations='Local optimization and one fixed architecture only. Components and velocity bounds differ from the 45-component architecture; this tests combined structure/initialization sensitivity. Covariance transfer and spectrograph calibration remain conditional.')
    (OUT/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--name');parser.add_argument('--free',action='store_true');parser.add_argument('--start');parser.add_argument('--max-nfev',type=int,default=400)
    args=parser.parse_args()
    if args.name:fit(args.name,args.free,args.start,args.max_nfev)
    else:workflow()
