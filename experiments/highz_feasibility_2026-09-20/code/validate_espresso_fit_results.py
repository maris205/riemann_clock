#!/usr/bin/env python3
"""Reconstruct saved ESPRESSO fit outputs and audit conditional comparisons.

No nonlinear optimization is rerun. Reported projected errors, when present,
are local diagnostics only; they do not resolve nuisance bounds or model choice.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import hashlib
import json
from pathlib import Path
import numpy as np
import espresso_conventional_null as e

def main():
    records=[];checks=[]
    def check(name,passed,**details):
        checks.append(dict(name=name,passed=bool(passed),**details))
        print(('PASS ' if passed else 'FAIL ')+name,flush=True)
    for path in sorted(e.OUT.glob('*.json')):
        if path.stem.endswith('_checkpoint'):continue
        d=json.loads(path.read_text())
        if not all(k in d for k in ['parameters','configuration','labels','free_shifts']):continue
        if not path.with_suffix('.npz').is_file():continue
        model=e.Model(free_shifts=d['free_shifts'],**d['configuration'])
        p=np.asarray(d['parameters'])
        residual,J,profiles=model.evaluate(p)
        check(path.stem+' parameter labels',model.labels==d['labels'])
        check(path.stem+' ndata',model.ndata==d['ndata'])
        delta=float(residual@residual-d['chi2'])
        check(path.stem+' objective reconstruction',abs(delta)<1e-8,
              chi2_difference=delta)
        point,lower,upper=model.initial()
        check(path.stem+' parameters in bounds',np.all(p>=lower-1e-9) and np.all(p<=upper+1e-9))
        saved=np.load(path.with_suffix('.npz'))
        check(path.stem+' saved residuals',np.allclose(saved['residuals'],residual,rtol=0,atol=1e-10))
        check(path.stem+' saved Jacobian',np.allclose(saved['jacobian'],J,rtol=1e-11,atol=1e-10))
        rank_values=np.linalg.svd(J,compute_uv=False)
        check(path.stem+' reported rank',int(np.sum(rank_values>rank_values[0]*1e-10))==d['jacobian_rank'])
        for L,profile,sl,line_d in zip(model.lines,profiles,model.slices,d['per_line']):
            key=str(L['key'])
            check(path.stem+' saved profile '+key,np.allclose(saved[key+'_model'],profile,rtol=0,atol=1e-12))
            check(path.stem+' line objective '+key,abs(residual[sl]@residual[sl]-line_d['chi2'])<1e-8)
            check(path.stem+' mask '+key,np.array_equal(saved[key+'_good'],L['good']))
        rec=dict(name=path.stem,optimizer_success=d['optimizer_success'],
                 optimizer_message=d['message'],optimality=d['optimality'],
                 nfev=d['nfev'],active_bounds=d['active_bounds'],chi2=d['chi2'],
                 chi2_per_nominal_ndf=d['chi2_per_nominal_ndf'],
                 condition_number=float(rank_values[0]/rank_values[-1]),
                 source_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        if model.free_shifts:
            shift_indices=[i for i,label in enumerate(model.labels) if label.startswith('shift_')]
            nuisance_indices=[i for i,label in enumerate(model.labels) if not label.startswith('shift_')]
            nuisance=J[:,nuisance_indices]
            shifts=J[:,shift_indices]
            U,sv,Vh=np.linalg.svd(nuisance,full_matrices=False)
            kept=sv>sv[0]*1e-10
            projected=shifts-U[:,kept]@(U[:,kept].T@shifts)
            ps=np.linalg.svd(projected,compute_uv=False)
            check(path.stem+' projected shift rank',np.sum(ps>ps[0]*1e-10)==len(model.shift_keys))
            covariance=np.linalg.inv(projected.T@projected)
            score=projected.T@residual
            remaining_shift=-covariance@score
            rec.update(projected_shift_singular_values=ps.tolist(),
                       shift_keys=model.shift_keys,
                       local_unbounded_shift_sigma_m_s=(1000*np.sqrt(np.diag(covariance))).tolist(),
                       local_unbounded_stationarity_shift_m_s=(1000*remaining_shift).tolist(),
                       local_unbounded_stationarity_delta_chi2=float(score@covariance@score),
                       shift_km_s=d['shifts_km_s'],
                       caution='Projected Gaussian errors ignore active bounds and model uncertainty; not discovery errors.')
        records.append(rec)
    out=dict(status='PASS' if checks and all(c['passed'] for c in checks) else 'NO_RESULTS' if not checks else 'FAIL',
             result_count=len(records),passed=sum(c['passed'] for c in checks),count=len(checks),
             checks=checks,results=records,
             model_source_sha256=hashlib.sha256(Path(e.__file__).read_bytes()).hexdigest(),
             validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             limitation='Saved-point replay and local rank audit, not independent nonlinear optimization or calibrated significance.')
    target=e.ROOT/'results/espresso_fit_results_validation.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['checks','results']},indent=2))
    if out['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
