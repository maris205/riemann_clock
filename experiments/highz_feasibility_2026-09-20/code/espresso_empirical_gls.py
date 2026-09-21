#!/usr/bin/env python3
"""Nonlinear re-fit under a covariance shape estimated outside absorption.

The continuum ACF is tapered and transferred to normalized absorption errors
as an explicit conditional assumption. The continuum variance rescaling is
reported separately; no calibrated physical-discovery probability is claimed.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import json,hashlib,time
import numpy as np
from scipy.linalg import cholesky,solve_triangular
from scipy.optimize import least_squares
from espresso_conventional_null import Model,ROOT
OUT=ROOT/'results/empirical_gls';BASE=ROOT/'results/espresso_null'
CONTROL=ROOT/'results/noise_covariance/controls.json'

class GLSModel(Model):
 def __init__(self,**kwargs):
  super().__init__(**kwargs)
  assert self.rho==0,'Do not double-whiten errors'
  controls=json.loads(CONTROL.read_text());acf=np.array(controls['primary']['acf'])
  self.kernel=acf*(1-np.arange(len(acf))/len(acf));self.variance=controls['primary']['variance']
  self.cholesky=[];self.white_cache=None
  for L in self.lines:
   idx=L['idx'][L['good']];lag=np.abs(idx[:,None]-idx[None,:])
   cov=np.zeros(lag.shape);valid=lag<len(self.kernel);cov[valid]=self.kernel[lag[valid]]
   self.cholesky.append(cholesky(cov,lower=True,check_finite=False))
 def evaluate(self,p,jac=True):
  if self.white_cache is not None and np.array_equal(p,self.white_cache[0]):return self.white_cache[1:]
  r,J,profiles=super().evaluate(p,jac)
  rw=np.empty_like(r);Jw=np.empty_like(J)
  for sl,chol in zip(self.slices,self.cholesky):
   rw[sl]=solve_triangular(chol,r[sl],lower=True,check_finite=False)
   Jw[sl]=solve_triangular(chol,J[sl],lower=True,check_finite=False)
  self.white_cache=(p.copy(),rw,Jw,profiles);return self.white_cache[1:]

def fit(name,free,start,max_nfev=350):
 OUT.mkdir(parents=True,exist_ok=True);m=GLSModel(free_shifts=free)
 p,lo,hi=m.initial();old=json.loads(Path(start).read_text());lookup=dict(zip(old['labels'],old['parameters']));p=np.array([lookup.get(k,v) for k,v in zip(m.labels,p)]);p=np.clip(p,lo+1e-9,hi-1e-9)
 initial=p.copy();began=time.time();calls=0
 def fun(p):
  nonlocal calls
  r=m.evaluate(p)[0];calls+=1
  if calls%25==0:
   print(name,calls,float(r@r),flush=True)
   (OUT/f'{name}_checkpoint.json').write_text(json.dumps(dict(labels=m.labels,parameters=p.tolist(),chi2=float(r@r))))
  return r
 print('START',name,m.ndata,m.npar,flush=True)
 opt=least_squares(fun,p,jac=m.jac,bounds=(lo,hi),x_scale='jac',ftol=1e-7,xtol=1e-8,gtol=1e-6,max_nfev=max_nfev)
 r,J,profiles=m.evaluate(opt.x);ss=np.linalg.svd(J,compute_uv=False);chi=float(r@r)
 result=dict(name=name,free_shifts=free,ndata=m.ndata,npar=m.npar,chi2=chi,nominal_ndf=m.ndata-m.npar,chi2_per_nominal_ndf=chi/(m.ndata-m.npar),optimizer_success=bool(opt.success),message=opt.message,nfev=opt.nfev,optimality=float(opt.optimality),elapsed_seconds=time.time()-began,labels=m.labels,parameters=opt.x.tolist(),initial_parameters=initial.tolist(),start_file=str(Path(start).relative_to(ROOT)),jacobian_rank=int(sum(ss>ss[0]*1e-10)),active_bounds=[m.labels[i] for i in range(m.npar) if min(opt.x[i]-lo[i],hi[i]-opt.x[i])<1e-5],shifts_m_s={k:1000*float(opt.x[m.labels.index(f'shift_{k}')]) for k in m.shift_keys},covariance_kernel=m.kernel.tolist(),continuum_variance=m.variance,chi2_if_continuum_variance_transferred=chi/m.variance,covariance_source_sha256=hashlib.sha256(CONTROL.read_bytes()).hexdigest(),interpretation='Constrained nonlinear GLS under a fixed continuum-informed stationary covariance shape. Covariance transfer and calibration/model uncertainty remain conditional.')
 (OUT/f'{name}.json').write_text(json.dumps(result,indent=2)+'\n');arr=dict(parameters=opt.x,residuals=r,jacobian=J)
 for L,profile in zip(m.lines,profiles):
  for k in ['v','wave','flux','error','good']:arr[f'{L["key"]}_{k}']=L[k]
  arr[f'{L["key"]}_model']=profile
 np.savez_compressed(OUT/f'{name}.npz',**arr);print('DONE',name,chi,opt.success,flush=True);return result

if __name__=='__main__':
 a=fit('tapered_null',False,BASE/'null_cross.json')
 if not a['optimizer_success']:a=fit('tapered_null_continued',False,OUT/(a['name']+'.json'))
 b=fit('tapered_alternative',True,OUT/(a['name']+'.json'))
 if not b['optimizer_success']:b=fit('tapered_alternative_continued',True,OUT/(b['name']+'.json'))
 c=fit('tapered_alternative_cross',True,BASE/'alternative_from_published.json')
 bs=[r for r in [b,c] if r['optimizer_success']];best=min(bs or [b,c],key=lambda x:x['chi2'])
 d=fit('tapered_null_cross',False,OUT/(best['name']+'.json'))
 aa=[r for r in [a,d] if r['optimizer_success']];an=min(aa or [a,d],key=lambda x:x['chi2'])
 summary=dict(null=an['name'],alternative=best['name'],null_success=an['optimizer_success'],alternative_success=best['optimizer_success'],ndata=an['ndata'],added_parameters=best['npar']-an['npar'],delta_chi2=an['chi2']-best['chi2'],delta_chi2_if_continuum_variance_transferred=(an['chi2']-best['chi2'])/an['continuum_variance'],shifts_m_s=best['shifts_m_s'],covariance_kernel=an['covariance_kernel'],continuum_variance=an['continuum_variance'])
 (OUT/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
