#!/usr/bin/env python3
"""Independent audit of new UVES stationarity and LSF-bound controls.

No fits are launched; older products are protected by pre-recorded hashes.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import hashlib,json,tempfile
import numpy as np
import uves_completion_fit as fit
from cross_instrument_validate import independent_profile

ROOT=fit.ROOT
OUT=ROOT/'results/uves_completion'
DEST=OUT/'validation.json'
CHECKS=[]

def check(name,passed,**details):
    CHECKS.append(dict(name=name,passed=bool(passed),**details))
    print(('PASS ' if passed else 'FAIL ')+name,flush=True)

def independent_stationarity(p,lower,upper,r,j):
    width=upper-lower
    gradient=np.array([np.dot(j[:,k],r) for k in range(len(p))])
    q=(p-lower)/width
    scaled_gradient=width*gradient
    mapping=np.empty(len(p))
    for k in range(len(p)):
        trial=q[k]-scaled_gradient[k]
        projection=min(1.,max(0.,trial))
        mapping[k]=q[k]-projection
    # Coleman--Li complementarity for TRF, in physical coordinates.
    complementarity=np.maximum((p-lower)*np.maximum(gradient,0.),
                               (upper-p)*np.maximum(-gradient,0.))
    return dict(mapping=float(np.max(np.abs(mapping))),
                coleman_li=float(np.max(complementarity)),
                normalized_primal_violation=float(max(np.max(-q),np.max(q-1.),0.)),
                worst_index=int(np.argmax(np.abs(mapping))),gradient=gradient,
                mapping_vector=mapping)

def main():
    previous=json.loads(DEST.read_text()) if DEST.exists() else {}
    protected=previous.get('protected_previous_files_initial_sha256',{})
    check('previous-file preservation baseline exists',len(protected)>=80,count=len(protected))
    for relative,wanted in protected.items():
        path=ROOT/relative
        current=hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
        check('preserved '+relative,current==wanted)
    # Exercise the actual stationarity function on hand-verifiable KKT cases.
    fixtures=[('interior optimum',.4,0.,0.),('lower correct',0.,2.,0.),
              ('upper correct',1.,-2.,0.),('lower wrong',0.,-2.,1.),
              ('upper wrong',1.,2.,1.)]
    for name,x,g,wanted in fixtures:
        diag=fit.stationarity(np.array([x]),np.array([0.]),np.array([1.]),np.array([g]),np.ones((1,1)))
        check('KKT fixture '+name,abs(diag['unit_box_gradient_mapping_inf']-wanted)<1e-14)
    p=np.array([1e-12]);r=np.array([1e12]);j=np.ones((1,1))
    edge=independent_stationarity(p,np.array([0.]),np.array([1.]),r,j)
    diag=fit.stationarity(p,np.array([0.]),np.array([1.]),r,j)
    check('near-bound mapping needs complementarity safeguard',
          diag['unit_box_gradient_mapping_inf']<1e-5 and edge['coleman_li']>1e-6,
          mapping=diag['unit_box_gradient_mapping_inf'],coleman_li=edge['coleman_li'])
    audited=[]
    with tempfile.TemporaryDirectory(prefix='uves_completion_audit_') as temporary:
        fit.OUT=Path(temporary)
        original=fit.make_model(False,False)
        wider=fit.make_model(False,True)
        op,ol,oh=original.initial();wp,wl,wh=wider.initial()
        n=original.base_npar
        check('wider bounds preserve labels and parameter count',original.labels==wider.labels and original.npar==147)
        check('wider bounds preserve every gas/continuum/zero bound',np.array_equal(ol[:n],wl[:n]) and np.array_equal(oh[:n],wh[:n]))
        check('only three appended LSF bounds widen',
              np.allclose(wl[n:],np.log(.5)) and np.allclose(wh[n:],np.log(1.4)) and
              np.allclose(ol[n:],np.log(.65)) and np.allclose(oh[n:],np.log(1.25)) and len(wl[n:])==3)
        check('wider bounds preserve physical residuals at shared parameters',np.array_equal(original.evaluate(op)[0],wider.evaluate(op)[0]))
        for wide in [False,True]:
            m0=fit.make_model(False,wide);m1=fit.make_model(True,wide)
            p0,l0,h0=m0.initial();p1,l1,h1=m1.initial()
            check(f'nuisance parity wide={wide}',set(m0.labels)==set(m1.labels)-{'shift_2382','shift_2600'})
            check(f'nested zero-shift residuals wide={wide}',np.array_equal(m0.evaluate(p0)[0],m1.evaluate(p1)[0]))
        for path in sorted(OUT.glob('*.json')):
            if path.name.endswith(('_start.json','_checkpoint.json')) or path.name=='validation.json':continue
            data=json.loads(path.read_text())
            if 'stationarity' not in data or 'configuration' not in data or 'parameters' not in data:continue
            if not path.with_suffix('.npz').exists():continue
            model=fit.make_model(data['free_shifts'],data['wide_lsf'],data['configuration']['oversample'])
            p=np.array(data['parameters']);r,j,profiles=model.evaluate(p)
            _,lo,hi=model.initial();arrays=np.load(path.with_suffix('.npz'))
            prefix=path.stem
            check(prefix+' source hashes',all((ROOT/relative).exists() and
                  hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==wanted
                  for relative,wanted in data['source_hashes_at_import'].items()))
            check(prefix+' labels',model.labels==data['labels'])
            check(prefix+' bounds',np.array_equal(lo,data['lower_bounds']) and np.array_equal(hi,data['upper_bounds']) and
                  np.array_equal(lo,arrays['lower']) and np.array_equal(hi,arrays['upper']))
            check(prefix+' accounting',model.ndata==data['ndata']==511 and model.npar==data['npar'])
            check(prefix+' reported LSF scales',all(abs(data['lsf_scale'][str(line['key'])]-
                  np.exp(p[model.base_npar+i]))<1e-12 for i,line in enumerate(model.lines)))
            check(prefix+' reported relative shifts',all(abs(data['shifts_km_s'][str(key)]-
                  p[model.labels.index(f'shift_{key}')])<1e-12 for key in model.shift_keys))
            check(prefix+' chi2',abs(float(r@r)-data['chi2'])<1e-7,recomputed=float(r@r),saved=data['chi2'])
            check(prefix+' residual and Jacobian arrays',np.max(np.abs(r-arrays['residuals']))<1e-8 and np.max(np.abs(j-arrays['jacobian']))<1e-7)
            computed=independent_stationarity(p,lo,hi,r,j)
            worst=computed['worst_index'];worst_label=model.labels[worst]
            step=5e-5 if worst_label.startswith(('v_','shift_')) else 1e-5
            plus=p.copy();minus=p.copy();plus[worst]+=step;minus[worst]-=step
            finite=(model.evaluate(plus)[0]-model.evaluate(minus)[0])/(2*step)
            derivative_error=float(np.linalg.norm(finite-j[:,worst])/max(np.linalg.norm(j[:,worst]),1e-14))
            check(prefix+' fitted worst-mapping derivative',derivative_error<7e-5,
                  parameter=worst_label,step=step,relative_l2_error=derivative_error,
                  finite_derivative_dot_residual=float(finite@r),analytic_gradient=float(computed['gradient'][worst]))
            # Restore per-line LSF state before independent physical profiles.
            model.evaluate(p)
            mapping=computed['mapping'];saved_mapping=data['stationarity']['unit_box_gradient_mapping_inf']
            check(prefix+' independent unit-box gradient mapping',abs(mapping-saved_mapping)<5e-7,
                  recomputed=mapping,saved=saved_mapping)
            check(prefix+' primal feasibility',computed['normalized_primal_violation']<1e-12,
                  normalized_violation=computed['normalized_primal_violation'])
            if data['termination_settings']['method']=='trf':
                recomputed=computed['coleman_li']
                check(prefix+' independent native TRF optimality',abs(recomputed-data['optimality'])<max(2e-7,abs(recomputed)*1e-5),
                      recomputed=recomputed,saved=data['optimality'])
            else:
                # Dogbox excludes exactly bound-active compatible coordinates.
                g=computed['gradient'].copy()
                g[((p==lo)&(g>0))|((p==hi)&(g<0))]=0.
                recomputed=float(np.max(np.abs(g)))
                check(prefix+' independent native dogbox optimality',abs(recomputed-data['optimality'])<max(2e-7,abs(recomputed)*1e-5),
                      recomputed=recomputed,saved=data['optimality'])
            classified=(data['optimality']<=1e-6 and mapping<=1e-5)
            check(prefix+' stationarity classification',classified==data['stationary_by_both_criteria'])
            snapshot=ROOT/data['initial_snapshot'];raw=snapshot.read_bytes();snap=json.loads(raw)
            check(prefix+' captured start snapshot hash',hashlib.sha256(raw).hexdigest()==data['initial_snapshot_sha256'])
            original_input=Path(snap['input_path'])
            check(prefix+' initial source hash',original_input.exists() and hashlib.sha256(original_input.read_bytes()).hexdigest()==snap['input_sha256'])
            for ki,(line,profile) in enumerate(zip(model.lines,profiles)):
                key=line['key'];independent=independent_profile(model,p,line,ki)
                error=float(np.max(np.abs(profile-independent)))
                check(prefix+f' independent fitted CGS profile {key}',error<2e-7,max_flux_error=error)
                check(prefix+f' saved pixel arrays {key}',all(np.array_equal(line[k],arrays[f'{key}_{k}']) for k in ['wave','v','flux','error','good']) and np.max(np.abs(profile-arrays[f'{key}_model']))<1e-10)
            fine=fit.make_model(data['free_shifts'],data['wide_lsf'],49)
            rf,_,_=fine.evaluate(p);delta=rf-r
            check(prefix+' 21-to-49 quadrature',float(delta@delta)<.05,
                  squared_noise_weighted_difference=float(delta@delta),max_sigma_difference=float(np.max(np.abs(delta))),
                  chi2_difference=float(rf@rf-r@r))
            singular=np.linalg.svd(j/np.maximum(np.linalg.norm(j,axis=0),1e-30),compute_uv=False)
            normalized_condition=float(singular[0]/singular[-1]) if singular[-1]>0 else None
            original_feasible=None;shared_chi2=None
            if data['wide_lsf'] and not data['free_shifts']:
                shared=fit.make_model(False,False,data['configuration']['oversample'])
                _,shared_lo,shared_hi=shared.initial()
                original_feasible=bool(np.all(p>=shared_lo) and np.all(p<=shared_hi))
                if original_feasible:
                    shared_r=shared.evaluate(p)[0];shared_chi2=float(shared_r@shared_r)
                    check(prefix+' feasible original-domain identical objective',
                          shared.labels==model.labels and np.array_equal(shared_r,r),
                          original_domain_chi2=shared_chi2,wide_domain_chi2=float(r@r))
            audited.append(dict(name=prefix,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),chi2=data['chi2'],
                optimizer_success=data['optimizer_success'],message=data['message'],native_optimality=data['optimality'],
                independently_computed_native_optimality=recomputed,unit_box_mapping=mapping,
                stationary_by_both_criteria=data['stationary_by_both_criteria'],
                most_nonstationary_parameter=data['labels'][computed['worst_index']],
                active_bounds=data['active_bounds'],lsf_scale=data['lsf_scale'],shifts_km_s=data['shifts_km_s'],
                column_normalized_jacobian_condition=normalized_condition,
                feasible_in_original_lsf_domain=original_feasible,original_domain_chi2=shared_chi2,
                free_shifts=data['free_shifts'],wide_lsf=data['wide_lsf']))
    missing=sorted({(False,False),(True,False),(False,True),(True,True)}-
                   {(a['free_shifts'],a['wide_lsf']) for a in audited})
    status=('FAIL' if not all(c['passed'] for c in CHECKS) else 'PENDING_NEW_FIT_PRODUCTS' if missing else 'PASS')
    out=dict(status=status,
        passed=sum(c['passed'] for c in CHECKS),count=len(CHECKS),checks=CHECKS,audited_results=audited,
        missing_configuration_pairs=[dict(free_shifts=free,wide_lsf=wide) for free,wide in missing],
        protected_previous_files_initial_sha256=protected,
        source_sha256=hashlib.sha256(Path(fit.__file__).read_bytes()).hexdigest(),
        validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='Verification of objective, bounds, stationarity diagnostics and saved arrays. PASS does not mean a fit meets stationarity thresholds or that a global/physical signal is established.')
    DEST.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','passed','count']},indent=2))
    if out['status']=='FAIL':raise SystemExit(1)

if __name__=='__main__':main()
