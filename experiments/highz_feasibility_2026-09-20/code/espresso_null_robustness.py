#!/usr/bin/env python3
"""Paired physical null/shift alternatives under predeclared stress scenarios."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from espresso_conventional_null import fit_case,OUT
CASES={'expected_fluctuation':dict(error_row=2),'ar1_rho03':dict(rho=.3),'lsf_minus3pct':dict(fwhm_scale=.97),'lsf_plus3pct':dict(fwhm_scale=1.03),'omit_blended':dict(exclude=[2344,2586])}
def run_case(item):
 name,kw=item
 # These initial robustness runs used the published component solution,
 # independently of the primary fit (preserve that initialization on replay).
 a=fit_case(name+'_null',start=None,max_nfev=350,**kw)
 b=fit_case(name+'_alternative',free=True,start=str(OUT/(name+'_null.json')),max_nfev=350,**kw)
 return name,a['chi2']-b['chi2']
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=5) as pool:
  for out in pool.map(run_case,CASES.items()):print('PAIRED_RESULT',out,flush=True)
