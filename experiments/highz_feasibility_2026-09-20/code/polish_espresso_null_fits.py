#!/usr/bin/env python3
"""Continue saved fit attempts that exhausted their evaluation budget.
Original attempts remain archived; the summary prefers converged saved results.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from concurrent.futures import ProcessPoolExecutor
import json
from espresso_conventional_null import fit_case,OUT
NAMES=['expected_fluctuation_null','omit_blended_alternative','omit_saturated_null']
def run(n):
 p=OUT/(n+'.json');x=json.loads(p.read_text())
 if x['optimizer_success']:return n,'already converged'
 out=fit_case(n+'_polished',free=x['free_shifts'],start=str(p),max_nfev=350,**x['configuration'])
 return n,out['optimizer_success'],out['chi2']
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=3) as pool:
  for r in pool.map(run,NAMES):print(r,flush=True)
