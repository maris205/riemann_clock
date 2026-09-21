#!/usr/bin/env python3
"""Exploratory linear opacity/shift sensitivity; not a nonlinear physical fit."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json
import numpy as np
from scipy.optimize import lsq_linear
import espresso_conventional_null as e

def main():
    base=json.loads((e.OUT/'null.json').read_text())
    model=e.Model(free_shifts=True,**base['configuration'])
    p,_,_=model.initial()
    lookup=dict(zip(base['labels'],base['parameters']))
    p=np.array([lookup.get(k,v) for k,v in zip(model.labels,p)])
    residual,J,_=model.evaluate(p)
    nuisance=J[:,:model.ngas+model.ncont]
    shifts=J[:,model.ngas+model.ncont:]
    U,s,Vh=np.linalg.svd(nuisance,full_matrices=False)
    U=U[:,s>s[0]*1e-10]
    project=lambda x:x-U@(U.T@x)
    rp=project(residual);Sp=project(shifts)
    opacity=np.zeros_like(shifts)
    # dF/dln(f) equals the sum of dF/dlog10(N_i), divided by ln(10),
    # within a transition. The reference line remains fixed to remove the
    # exact common-opacity/column-density degeneracy.
    for key,sl in zip(model.keys,model.slices):
        if key in model.shift_keys:
            opacity[sl,model.shift_keys.index(key)]=J[sl,:model.ncomp].sum(axis=1)/np.log(10)
    Ap=project(opacity)
    cases=[('shifts',Sp,-.5,.5),('opacity5pct',Ap,-.05,.05),
           ('opacity10pct',Ap,-.1,.1),
           ('shifts_and_opacity10pct',np.column_stack([Sp,Ap]),
            np.r_[np.full(5,-.5),np.full(5,-.1)],np.r_[np.full(5,.5),np.full(5,.1)])]
    results=[]
    for label,X,lower,upper in cases:
        fit=lsq_linear(X,-rp,bounds=(lower,upper),tol=1e-12)
        results.append(dict(label=label,local_projected_delta_chi2=float(
            rp@rp-np.linalg.norm(rp+X@fit.x)**2),parameters=fit.x.tolist()))
    out=dict(caution='Local linear sensitivity only; nuisance directions unbounded; not a nonlinear fit or measured laboratory prior.',
             shift_keys=model.shift_keys,results=results)
    (e.ROOT/'results/espresso_opacity_local_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
