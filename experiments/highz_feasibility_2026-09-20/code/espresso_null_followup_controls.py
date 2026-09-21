#!/usr/bin/env python3
"""Post-inspection conventional controls motivated by strong-line shift pattern.
These are exploratory stress ranges, not externally measured calibration priors.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from concurrent.futures import ProcessPoolExecutor
from espresso_conventional_null import fit_case,OUT
CASES={'opacity_zero':dict(strength_nuisance=True,zero_nuisance=True),'omit_saturated':dict(exclude=[2382,2600]),'wider_bounds':dict(velocity_span=5.,continuum_span=.05)}
def run_case(item):
 name,kw=item
 a=fit_case(name+'_null',start=str(OUT/'null.json'),max_nfev=350,**kw)
 b=fit_case(name+'_alternative',free=True,start=str(OUT/(name+'_null.json')),max_nfev=350,**kw)
 return name,a['chi2']-b['chi2']
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=3) as pool:
  for out in pool.map(run_case,CASES.items()):print('PAIRED_RESULT',out,flush=True)
