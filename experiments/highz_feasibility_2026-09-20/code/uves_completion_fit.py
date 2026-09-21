#!/usr/bin/env python3
"""Bounded UVES endpoint completion and paired wider Gaussian-LSF sensitivity.

Preserves earlier files. Expected-fluctuation errors, three clean FeII windows,
same gas/continuum/zero nuisance under both hypotheses. A successful ftol/xtol
exit is not labelled first-order stationarity. No discovery significance.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import argparse,datetime,hashlib,json,time
import numpy as np
from scipy.optimize import least_squares
import scipy
import cross_instrument_uves as uv
ROOT=uv.ROOT;OUT=ROOT/'results/uves_completion'
CODE_HASH=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
SOURCES={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'code/cross_instrument_uves.py',ROOT/'code/espresso_conventional_null.py']}

class CompletionModel(uv.FlexibleLSFModel):
    def __init__(self,*args,lsf_bounds=(.65,1.25),**kw):
        self.lsf_bounds=tuple(lsf_bounds)
        super().__init__(*args,**kw)
    def initial(self):
        p,lo,hi=super().initial()
        lo[self.base_npar:]=np.log(self.lsf_bounds[0]);hi[self.base_npar:]=np.log(self.lsf_bounds[1])
        return p,lo,hi

def make_model(free=False,wide=False,oversample=21):
    OUT.mkdir(parents=True,exist_ok=True)
    uv.OUT=OUT  # isolate load_uves metadata side effects from earlier results
    uv.configure()
    return CompletionModel(free_shifts=free,continuum=True,error_row=2,oversample=oversample,
       exclude=[2260,2344,2586],zero_nuisance=True,velocity_span=5.,continuum_span=.03,
       lsf_bounds=(.5,1.4) if wide else (.65,1.25))

def stationarity(p,lo,hi,r,J):
    # Box-normalized gradient mapping for f = (1/2) ||r||^2, independent of SciPy.
    span=hi-lo;q=(p-lo)/span;g=J.T@r;gq=span*g
    mapping=q-np.clip(q-gq,0,1)
    near_low=(p-lo)<=1e-7*np.maximum(1.,span)
    near_high=(hi-p)<=1e-7*np.maximum(1.,span)
    projected=g.copy();projected[near_low&(g>=0)]=0;projected[near_high&(g<=0)]=0
    norm=np.linalg.norm(J,axis=0)
    return dict(objective='0.5 * sum squared weighted residuals',
       unit_box_gradient_mapping_inf=float(np.max(abs(mapping))),
       unit_box_mapping_definition='q=(p-lo)/(hi-lo); gq=(hi-lo)*(J.T@r); infinity norm of q-clip(q-gq,0,1)',
       raw_gradient_inf=float(np.max(abs(g))),projected_gradient_inf=float(np.max(abs(projected))),
       projected_column_normalized_gradient_inf=float(np.max(abs(projected)/np.maximum(norm,1e-14))),
       max_bound_violation=float(max(np.max(lo-p),np.max(p-hi),0)),
       nearest_normalized_bound_distance=float(np.min(np.minimum(q,1-q))),
       near_lower=np.flatnonzero(near_low).tolist(),near_upper=np.flatnonzero(near_high).tolist(),
       most_nonstationary_parameter_index=int(np.argmax(abs(mapping))))

def run(name,start,free=False,wide=False,max_nfev=1800,oversample=21,xscale='jac',method='trf'):
    OUT.mkdir(parents=True,exist_ok=True)
    m=make_model(free,wide,oversample);p,lo,hi=m.initial()
    raw=Path(start).read_bytes();old=json.loads(raw);lookup=dict(zip(old['labels'],old['parameters']))
    p=np.array([lookup.get(k,val) for k,val in zip(m.labels,p)]);p=np.clip(p,lo+1e-10,hi-1e-10)
    snapshot=dict(input_path=str(start),input_sha256=hashlib.sha256(raw).hexdigest(),
        original_record=old,labels=m.labels,parameters=p.tolist(),lower=lo.tolist(),upper=hi.tolist(),
        captured_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    snapshot_path=OUT/(name+'_start.json');snapshot_path.write_text(json.dumps(snapshot,indent=2)+'\n')
    t=time.time();counter=0;best=(np.inf,p.copy());trace=[]
    def fun(x):
        nonlocal counter,best
        r=m.evaluate(x)[0];chi=float(r@r);counter+=1
        if chi<best[0]:best=(chi,x.copy())
        if counter%25==0:
            state=dict(labels=m.labels,parameters=best[1].tolist(),chi2=best[0],evaluation=counter,
                last_trial_chi2=chi,seconds=time.time()-t)
            (OUT/(name+'_checkpoint.json')).write_text(json.dumps(state,indent=2)+'\n')
            trace.append({k:state[k] for k in ['evaluation','chi2','last_trial_chi2','seconds']})
            print(f'{name}: nfev {counter}, best chi2 {best[0]:.9f}, trial {chi:.9f}',flush=True)
        return r
    r0,J0,_=m.evaluate(p);initial_stat=stationarity(p,lo,hi,r0,J0)
    print(f'{name}: ndata={m.ndata}, npar={m.npar}, initial chi2={r0@r0:.9f}, mapping={initial_stat["unit_box_gradient_mapping_inf"]:.6g}',flush=True)
    if xscale=='physical':
        scale=np.array([.2 if k.startswith(('logN_','v_','logb_')) else .003 if k.startswith('continuum') else .002 if k.startswith('zero_') else .05 for k in m.labels])
    else:scale='jac'
    settings=dict(method=method,ftol=1e-10,xtol=1e-11,gtol=1e-6,max_nfev=max_nfev,x_scale=xscale)
    opt=least_squares(fun,p,jac=m.jac,bounds=(lo,hi),x_scale=scale,method=method,
        ftol=settings['ftol'],xtol=settings['xtol'],gtol=settings['gtol'],max_nfev=max_nfev)
    r,J,models=m.evaluate(opt.x);st=stationarity(opt.x,lo,hi,r,J)
    sv=np.linalg.svd(J,compute_uv=False)
    criteria=dict(scipy_scaled_gradient_below_1e_6=bool(opt.optimality<=1e-6),unit_box_mapping_below_1e_5=bool(st['unit_box_gradient_mapping_inf']<=1e-5))
    out=dict(name=name,free_shifts=free,wide_lsf=wide,lsf_bounds=list(m.lsf_bounds),
        configuration=dict(error_row=2,oversample=oversample,exclude=[2260,2344,2586],zero_nuisance=True,velocity_span=5.,continuum_span=.03),
        ndata=m.ndata,npar=m.npar,chi2=float(r@r),labels=m.labels,parameters=opt.x.tolist(),lower_bounds=lo.tolist(),upper_bounds=hi.tolist(),
        optimizer_success=bool(opt.success),optimizer_status=int(opt.status),message=str(opt.message),
        nfev=int(opt.nfev),njev=int(opt.njev),optimality=float(opt.optimality),seconds=time.time()-t,
        termination_settings=settings,stationarity=st,stationarity_criteria=criteria,stationary_by_both_criteria=bool(all(criteria.values())),
        numerical_rank=int(sum(sv>sv[0]*1e-10)),singular_values=sv.tolist(),
        active_bounds=[m.labels[i] for i in range(m.npar) if min(opt.x[i]-lo[i],hi[i]-opt.x[i])<=1e-7*max(1,hi[i]-lo[i])],
        shifts_km_s={k:float(opt.x[m.labels.index(f'shift_{k}')]) for k in m.shift_keys},
        lsf_scale={L['key']:float(np.exp(opt.x[m.base_npar+i])) for i,L in enumerate(m.lines)},
        initial_chi2=float(r0@r0),initial_stationarity=initial_stat,initial_snapshot=str(snapshot_path.relative_to(ROOT)),
        initial_snapshot_sha256=hashlib.sha256(snapshot_path.read_bytes()).hexdigest(),source_hashes_at_import=SOURCES,
        versions=dict(numpy=np.__version__,scipy=scipy.__version__),best_trial_chi2=best[0],trace=trace,
        scope='Conditional UVES3line joint gas test; no calibrated global significance and no cosmic-time law.')
    (OUT/(name+'.json')).write_text(json.dumps(out,indent=2)+'\n')
    arrays=dict(parameters=opt.x,lower=lo,upper=hi,residuals=r,jacobian=J)
    for L,mod in zip(m.lines,models):
        for key in ['v','wave','flux','error','good']:arrays[f'{L["key"]}_{key}']=L[key]
        arrays[f'{L["key"]}_model']=mod
    np.savez_compressed(OUT/(name+'.npz'),**arrays)
    print(json.dumps({k:out[k] for k in ['name','chi2','nfev','optimizer_success','optimizer_status','optimality','stationarity','stationarity_criteria','shifts_km_s','lsf_scale','active_bounds']},indent=2),flush=True)
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--name',required=True);p.add_argument('--start',required=True);p.add_argument('--free',action='store_true');p.add_argument('--wide',action='store_true');p.add_argument('--max-nfev',type=int,default=1800);p.add_argument('--oversample',type=int,default=21);p.add_argument('--xscale',choices=['jac','physical'],default='jac');p.add_argument('--method',choices=['trf','dogbox'],default='trf')
    a=p.parse_args();run(a.name,a.start,a.free,a.wide,a.max_nfev,a.oversample,a.xscale,a.method)
