#!/usr/bin/env python3
"""Frozen model-defined strong-core masks, with both hypotheses re-fitted.

Exploratory follow-up to the relative line-offset pattern, not a blinded test.
Masks depend on a fixed earlier null spectrum, never the current noisy flux or
moving fit. Original data and previously delivered fits are not modified.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import json,time,hashlib
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import least_squares
from scipy.ndimage import binary_dilation
from espresso_conventional_null import Model,ROOT
BASE=ROOT/'results/espresso_null';OUT=ROOT/'results/saturation_controls'

class CoreMaskModel(Model):
 def __init__(self,threshold=.1,padding_km_s=0.,**kwargs):
  super().__init__(**kwargs)
  frozen=np.load(BASE/'null_cross.npz');self.mask_audit=[];self.slices=[];n=0
  for L in self.lines:
   k=L['key'];original=L['good'].copy();remove=np.zeros(len(original),bool)
   assert np.allclose(L['wave'],frozen[f'{k}_wave'],rtol=0,atol=1e-9)
   if k in [2382,2600]:
    remove=frozen[f'{k}_model']<threshold
    dilation=int(np.ceil(padding_km_s/self.dv))
    if dilation:remove=binary_dilation(remove,iterations=dilation)
    L['good']=original & ~remove
   self.mask_audit.append(dict(line=k,threshold=threshold if k in [2382,2600] else None,padding_km_s=padding_km_s,original_pixels=int(original.sum()),removed_pixels=int((original&remove).sum()),retained_pixels=int(L['good'].sum())))
   self.slices.append(slice(n,n+int(L['good'].sum())));n+=int(L['good'].sum())
  self.ndata=n;self.cache=None

def fit(name,free,threshold,padding,start,max_nfev=350):
 OUT.mkdir(parents=True,exist_ok=True);m=CoreMaskModel(threshold=threshold,padding_km_s=padding,free_shifts=free)
 p,lo,hi=m.initial();old=json.loads(Path(start).read_text());lookup=dict(zip(old['labels'],old['parameters']));p=np.array([lookup.get(k,v) for k,v in zip(m.labels,p)])
 p=np.clip(p,lo+1e-9,hi-1e-9);initial=p.copy();began=time.time();calls=0
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
 result=dict(name=name,free_shifts=free,threshold=threshold,padding_km_s=padding,ndata=m.ndata,npar=m.npar,chi2=chi,nominal_ndf=m.ndata-m.npar,chi2_per_nominal_ndf=chi/(m.ndata-m.npar),optimizer_success=bool(opt.success),message=opt.message,nfev=opt.nfev,optimality=float(opt.optimality),elapsed_seconds=time.time()-began,labels=m.labels,parameters=opt.x.tolist(),initial_parameters=initial.tolist(),start_file=str(Path(start).relative_to(ROOT)),mask_audit=m.mask_audit,jacobian_rank=int(sum(ss>ss[0]*1e-10)),active_bounds=[m.labels[i] for i in range(m.npar) if min(opt.x[i]-lo[i],hi[i]-opt.x[i])<1e-5],shifts_m_s={k:1000*float(opt.x[m.labels.index(f'shift_{k}')]) for k in m.shift_keys},mask_source_file='results/espresso_null/null_cross.npz',mask_source_sha256=hashlib.sha256((BASE/'null_cross.npz').read_bytes()).hexdigest(),interpretation='Frozen model-defined mask; joint bounded nonlinear refit; no complete covariance or discovery calibration.')
 (OUT/f'{name}.json').write_text(json.dumps(result,indent=2)+'\n');arr=dict(parameters=opt.x,residuals=r,jacobian=J)
 for L,profile in zip(m.lines,profiles):
  for k in ['v','wave','flux','error','good']:arr[f'{L["key"]}_{k}']=L[k]
  arr[f'{L["key"]}_model']=profile
 np.savez_compressed(OUT/f'{name}.npz',**arr)
 print('DONE',name,chi,opt.success,flush=True);return result

def run_case(case):
 name,threshold,padding=case
 a=fit(name+'_null',False,threshold,padding,BASE/'null_cross.json')
 if not a['optimizer_success']:a=fit(name+'_null_continued',False,threshold,padding,OUT/f'{name}_null.json')
 b=fit(name+'_alternative',True,threshold,padding,OUT/(a['name']+'.json'))
 if not b['optimizer_success']:b=fit(name+'_alternative_continued',True,threshold,padding,OUT/f'{name}_alternative.json')
 # A shift fit initialized from the full-data alternative provides a second
 # local basin check on the same fixed mask.
 c=fit(name+'_alternative_cross',True,threshold,padding,BASE/'alternative_from_published.json')
 alternatives=[z for z in [b,c] if z['optimizer_success']];best=min(alternatives or [b,c],key=lambda z:z['chi2'])
 # Cross-start the null from the selected shift fit, with shifts removed.
 d=fit(name+'_null_cross',False,threshold,padding,OUT/(best['name']+'.json'))
 nulls=[z for z in [a,d] if z['optimizer_success']];an=min(nulls or [a,d],key=lambda z:z['chi2'])
 result=dict(case=name,null=an['name'],alternative=best['name'],ndata=an['ndata'],added_parameters=best['npar']-an['npar'],delta_chi2=an['chi2']-best['chi2'],null_success=an['optimizer_success'],alternative_success=best['optimizer_success'],shifts_m_s=best['shifts_m_s'],mask_audit=an['mask_audit'])
 (OUT/f'{name}_comparison.json').write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=2) as pool:
  for result in pool.map(run_case,[('core01',.1,0.),('core03_padded',.3,1.)]):print(json.dumps(result),flush=True)
