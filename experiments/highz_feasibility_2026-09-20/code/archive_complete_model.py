#!/usr/bin/env python3
"""Bounded completion campaign for four previously screened FeII absorbers.

New adapter over the audited archive Voigt evaluator; previous code/results are
read-only. Every fit uses full-coadd native pixels and expected-fluctuation
errors. Component and window choices are exploratory and recorded before any
relative-shift fit. No alpha parameter or time law is fitted.
"""
from __future__ import annotations
import argparse, concurrent.futures as cf, hashlib, json, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import svd
import archive_profile_pilot as old

ROOT=old.ROOT; DATA=old.DATA; OUT=ROOT/'results/archive_complete'; C=old.C
NUISANCE={'continuum_constant':.10,'continuum_slope_per_100km_s':.03,
          'zero':.02,'lsf_factors':[.6,1.4],'b_km_s':[.5,30.],
          'log10_column':[9.,17.],'shift_km_s':[-1.,1.]}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_json(path):return json.loads(Path(path).read_text())

class Model(old.Model):
 def __init__(self,configuration,ncomp,free=False,oversample=9):
  self.cfg=dict(configuration);self.target=self.cfg['target'];self.ncomp=ncomp;self.free=free;self.oversample=oversample
  self.keys=list(map(int,self.cfg['lines']));assert 2374 in self.keys and 2382 in self.keys
  self.shiftkeys=[k for k in self.keys if k!=2374] if free else []
  self.labels=[f'logN_{i}' for i in range(ncomp)]+[f'v_{i}' for i in range(ncomp)]+[f'logb_{i}' for i in range(ncomp)]
  self.labels += [f'{a}_{k}' for k in self.keys for a in ['cont0','cont1','zero','logfwhm']]
  self.labels += [f'shift_{k}' for k in self.shiftkeys]
  self.npar=len(self.labels);self.lines=[];self.slices=[];self.cache=None;off=0
  d=np.load(DATA/f'{self.target}.npz');meta=read_json(DATA/f'{self.target}_metadata.json');ad=old.atoms(self.keys)
  wavelength=d['wavelength_AA'];dv=C*np.log(10)*meta['header']['CD1_1']
  for k in self.keys:
   a=ad[k];assert len(a)>0
   ref=np.sum(a[:,0]*a[:,1])/a[:,1].sum();vfull=C*np.log(wavelength/(ref*(1+self.cfg['z'])))
   use=(vfull>=self.cfg['window'][0])&(vfull<=self.cfg['window'][1]);ix=np.flatnonzero(use);v=vfull[ix]
   if len(ix)<3:raise ValueError(f'No complete native window: {self.target} {k}')
   good=d['valid'][ix].copy();error=d['expected_fluctuation'][ix]
   if not np.all(np.isfinite(error[good])&(error[good]>0)):raise ValueError('Expected errors unavailable on fixed source-valid pixels')
   row=min(meta['tables']['1']['rows'],key=lambda r:abs(r['Wavelength']-np.median(wavelength[ix])))
   nom=C/row['NomResolPower'];pad=int(np.ceil(6*nom*NUISANCE['lsf_factors'][1]/dv))+2
   grid=v[0]-pad*dv+(np.arange((len(v)+2*pad)*oversample)+.5)*dv/oversample-.5*dv
   self.lines.append(dict(key=k,atom=a,ref=ref,v=v,dv=dv,grid=grid,pad=pad,good=good,
      flux=d['flux'][ix],error=error,statistical_error=d['error'][ix],wave=wavelength[ix],
      source_indices=ix,x=(v-v.mean())/100,nominal_fwhm=nom,
      nominal_R=float(row['NomResolPower']),arc_R=float(row['ArcResolPower'])))
   self.slices.append(slice(off,off+good.sum()));off+=good.sum()
  self.ndata=int(off)

 def initial(self,seed=0):
  n=self.ncomp;rng=np.random.default_rng(260921+seed)
  mapping=self.cfg.get('centers_by_n',self.cfg.get('centers_by_count',{}))
  if str(n) in mapping:v=np.asarray(mapping[str(n)],float)
  else:
   pool=np.array(self.cfg.get('centers',np.linspace(*self.cfg['component_bounds'],n)),float)
   v=np.interp(np.linspace(0,len(pool)-1,n),np.arange(len(pool)),np.sort(pool))
  if len(v)!=n:raise ValueError('Wrong count in frozen center map')
  if seed:v=v+rng.normal(0,2.5,n)
  anchor=next(L for L in self.lines if L['key']==2374)
  weight=.08+np.maximum(0,1-np.interp(v,anchor['v'],anchor['flux']))
  logN=self.cfg.get('logtotal',15.2)+np.log10(weight/weight.sum())+rng.normal(0,.10,n)
  b=np.full(n,5.) if seed!=2 else np.full(n,9.)
  p=list(logN)+list(v)+list(np.log(b));lo=[9.]*n+[self.cfg['component_bounds'][0]]*n+[np.log(.5)]*n
  hi=[17.]*n+[self.cfg['component_bounds'][1]]*n+[np.log(30.)]*n
  for L in self.lines:
   side=L['good']&(np.abs(L['x'])>=np.quantile(np.abs(L['x']),.75))
   c0=np.clip(np.median(L['flux'][side])-1,-.09,.09) if side.any() else 0
   p += [c0,0,0,np.log(L['nominal_fwhm'])]
   lo += [-.10,-.03,-.02,np.log(.6*L['nominal_fwhm'])]
   hi += [.10,.03,.02,np.log(1.4*L['nominal_fwhm'])]
  p += [0]*len(self.shiftkeys);lo += [-1.]*len(self.shiftkeys);hi += [1.]*len(self.shiftkeys)
  return np.asarray(p),np.asarray(lo),np.asarray(hi)

def stationarity(p,lo,hi,r,J,optimality):
 g=J.T@r;s=hi-lo;q=(p-lo)/s;gq=g*s
 mapping=q-np.clip(q-gq,0,1)
 cl=np.max(np.maximum(q*np.maximum(gq,0),(1-q)*np.maximum(-gq,0)))
 norms=np.linalg.norm(J,axis=0);scaled=g/np.maximum(norms,1e-30)
 return {'scipy_optimality':float(optimality),'independent_coleman_li':float(cl),
         'normalized_gradient_mapping_inf':float(np.max(np.abs(mapping))),
         'column_normalized_gradient_inf':float(np.max(np.abs(scaled))),
         'interpretation':'Recorded first-order diagnostics; ftol/xtol termination does not prove a global optimum or eliminate weak/bound-active nuisance directions.'}

def fit(config,ncomp,free,seed=0,start=None,name=None,oversample=9,max_nfev=350,fixed_shift=None):
 m=Model(config,ncomp,free,oversample);p,lo,hi=m.initial(seed)
 if start is not None:
  lookup=dict(zip(start['labels'],start['parameters']));p=np.asarray([lookup.get(k,v) for k,v in zip(m.labels,p)])
 p=np.clip(p,lo+1e-8,hi-1e-8);initial=p.copy();output=OUT/config['target'];output.mkdir(parents=True,exist_ok=True)
 name=name or f'n{ncomp}_{"shift" if free else "null"}_s{seed}_os{oversample}'
 counter=0
 # Optional fixed 2382 displacement is a profile objective, not a pseudo-observation.
 fixed_index=m.labels.index('shift_2382') if fixed_shift is not None else None
 active=np.ones(m.npar,bool)
 if fixed_index is not None:active[fixed_index]=False;p[fixed_index]=fixed_shift
 def expand(q):
  full=p.copy();full[active]=q;return full
 def fun(q):
  nonlocal counter
  full=expand(q);r=m.fun(full);counter+=1
  if counter%50==0:
   (output/f'{name}_checkpoint.json').write_text(json.dumps({'name':name,'evaluation':counter,'chi2':float(r@r),'labels':m.labels,'parameters':full.tolist()},indent=2)+'\n')
   print(config['target'],name,counter,float(r@r),flush=True)
  return r
 def jac(q):return m.jac(expand(q))[:,active]
 tick=time.time();opt=least_squares(fun,p[active],jac=jac,bounds=(lo[active],hi[active]),x_scale='jac',
     ftol=2e-8,xtol=1e-8,gtol=1e-6,max_nfev=max_nfev)
 best=expand(opt.x);r,J,profiles=m.evaluate(best);chi=float(r@r);k=int(active.sum());ndf=m.ndata-k
 record={'name':name,'target':config['target'],'configuration':config,'ncomp':ncomp,'free_shifts':free,
  'seed':seed,'oversample':oversample,'max_nfev':max_nfev,'ndata':m.ndata,'npar':k,'nominal_ndf':ndf,
  'chi2':chi,'chi2_per_ndf':chi/ndf,'aicc':chi+2*k+2*k*(k+1)/(m.ndata-k-1),
  'optimizer_success':bool(opt.success),'optimizer_status':int(opt.status),'message':opt.message,
  'nfev':int(opt.nfev),'seconds':time.time()-tick,'labels':m.labels,'parameters':best.tolist(),
  'initial_parameters':initial.tolist(),'bounds_lower':lo.tolist(),'bounds_upper':hi.tolist(),
  'active_bounds':[m.labels[i] for i in range(m.npar) if min(best[i]-lo[i],hi[i]-best[i])<1e-5],
  'shifts_m_s':{str(k):1000*float(best[m.labels.index(f'shift_{k}')]) for k in m.shiftkeys},
  'fixed_2382_shift_km_s':fixed_shift,'error_kind':'FITS primary zero-based row2 normalized expected fluctuation; diagonal weighting',
  'stationarity':stationarity(best[active],lo[active],hi[active],r,J[:,active],opt.optimality),
  'per_line':[{'line':L['key'],'ndata':int(L['good'].sum()),'chi2':float(r[sl]@r[sl]),'chi2_per_pixel':float(r[sl]@r[sl]/L['good'].sum())} for L,sl in zip(m.lines,m.slices)],
  'source_hashes':{'coadd':sha(DATA/f'{config["target"]}.npz'),'metadata':sha(DATA/f'{config["target"]}_metadata.json'),
      'atomic':sha(old.ATOMIC),'base_model':sha(old.__file__),'adapter':sha(__file__)}}
 (output/f'{name}.json').write_text(json.dumps(record,indent=2)+'\n')
 arrays={'parameters':best,'residual':r,'jacobian':J,'fitted_parameter_mask':active}
 for L,profile in zip(m.lines,profiles):
  for key in ['v','wave','flux','error','statistical_error','good','source_indices']:arrays[f'{L["key"]}_{key}']=L[key]
  arrays[f'{L["key"]}_model']=profile
 np.savez_compressed(output/f'{name}.npz',**arrays)
 print(config['target'],name,'DONE',chi,chi/ndf,opt.success,flush=True)
 return record

def _job(args):return fit(**args)

def local_pair_uncertainty(record):
 folder=OUT/record['target'];saved=np.load(folder/f'{record["name"]}.npz');J=saved['jacobian']
 j=record['labels'].index('shift_2382');a=J[:,j];N=np.delete(J,j,axis=1)
 scale=np.linalg.norm(N,axis=0);Ns=N/np.maximum(scale,1e-100);u,s,vh=svd(Ns,full_matrices=False,lapack_driver='gesvd')
 answers=[]
 for cutoff in [1e-8,1e-10,1e-12]:
  rank=int(np.sum(s>s[0]*cutoff));perp=a-u[:,:rank]@(u[:,:rank].T@a);information=float(perp@perp)
  answers.append({'svd_relative_cutoff':cutoff,'nuisance_rank':rank,'information_per_km_s_squared':information,
      'conditional_sigma_m_s':1000/np.sqrt(information) if information>1e-20 else None})
 return {'observable':'velocity(FeII2382)-velocity(FeII2374)','shift_m_s':record['shifts_m_s']['2382'],
   'method':'Weighted-Jacobian 2382 column projected off ALL other gas/continuum/zero/LSF/shift nuisance columns, after nuisance-column normalization.',
   'rank_sensitivity':answers,'active_bounds':record['active_bounds'],
   'limitations':'Local linear, fixed selected architecture and expected diagonal errors; ignores active-bound truncation, calibration, blends and architecture uncertainty. Not alpha or a cosmic-time observation.'}

def run_campaign(configfile,workers=3):
 config=read_json(configfile);folder=OUT/config['target'];folder.mkdir(parents=True,exist_ok=True)
 counts=config.get('component_counts',[4,6,8,10,12,16]);seeds=config.get('seeds',[0,1,2]);budget=config.get('max_nfev',350)
 protocol={'frozen_utc':datetime.now(timezone.utc).isoformat(),'configuration_file':str(Path(configfile).resolve()),
   'configuration_sha256':sha(configfile),'nuisance_bounds':NUISANCE,'counts':counts,'seeds':seeds,'initial_max_nfev':budget,
   'initial_oversample':9,'final_oversample':21,'selection':'Minimum null AICc among saved successful stops. Also retain best incomplete endpoint. No component selection by displacement or free-shift improvement.',
   'adequacy_gate':'Conditional diagnostic only: null chi2/nominaldof<=1.5 and every line chi2/pixel<=1.8, successful optimizer termination. Does not prove model validity.',
   'bounded_continuation':'One extra null continuation for each best count if initial best stop unsuccessful; chosen-null OS21; H1 from that null and its seed; H0 cross from H1. Fixed-profile at most7points for2382 if possible.'}
 (folder/'campaign_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
 records=[]
 jobs=[dict(config=config,ncomp=n,free=False,seed=s,oversample=9,max_nfev=budget) for n in counts for s in seeds]
 if workers>1:
  with cf.ProcessPoolExecutor(max_workers=workers) as pool:
   for r in pool.map(_job,jobs):records.append(r)
 else:
  records=[_job(j) for j in jobs]
 bycount=[]
 for n in counts:
  opts=[r for r in records if r['ncomp']==n];best=min(opts,key=lambda r:r['chi2'])
  if not best['optimizer_success']:
   continued=fit(config,n,False,best['seed'],best,name=f'n{n}_null_continued_os9',oversample=9,max_nfev=budget)
   records.append(continued);opts.append(continued)
  good=[r for r in opts if r['optimizer_success']];chosen=min(good or opts,key=lambda r:r['chi2']);bycount.append(chosen)
 completed=[r for r in bycount if r['optimizer_success']]
 selected=min(completed or bycount,key=lambda r:r['aicc']);n=selected['ncomp']
 refined=fit(config,n,False,selected['seed'],selected,name=f'n{n}_null_refined_os21',oversample=21,max_nfev=500)
 ordinary={'target':config['target'],'by_count':[{'ncomp':r['ncomp'],'name':r['name'],'chi2':r['chi2'],'nominal_ndf':r['nominal_ndf'],'aicc':r['aicc'],'optimizer_success':r['optimizer_success'],'per_line':r['per_line']} for r in bycount],
    'selected_null':refined['name'],'ncomp':n,'selection_reason':protocol['selection'],
    'best_incomplete_endpoint':min([r for r in records if not r['optimizer_success']],key=lambda r:r['chi2'])['name'] if any(not r['optimizer_success'] for r in records) else None,
    'passes_conditional_adequacy_gate':bool(refined['optimizer_success'] and refined['chi2_per_ndf']<=1.5 and max(x['chi2_per_pixel'] for x in refined['per_line'])<=1.8)}
 (folder/'ordinary_model_summary.json').write_text(json.dumps(ordinary,indent=2)+'\n')
 if not ordinary['passes_conditional_adequacy_gate']:
  print(config['target'],'ordinary model inadequate or optimization incomplete; alternatives not auto-run',flush=True)
  return ordinary
 return run_alternative(config,refined,ordinary)

def run_alternative(config,refined,ordinary):
 folder=OUT/config['target'];n=refined['ncomp']
 alt=fit(config,n,True,refined['seed'],refined,name=f'n{n}_shift_fromnull_os21',oversample=21,max_nfev=600)
 cross=fit(config,n,False,refined['seed'],alt,name=f'n{n}_null_cross_os21',oversample=21,max_nfev=500)
 null=min([refined,cross],key=lambda r:r['chi2'])
 if null['chi2']<refined['chi2']-1e-3:
  alt2=fit(config,n,True,null['seed'],null,name=f'n{n}_shift_cross_os21',oversample=21,max_nfev=500)
  alt=min([alt,alt2],key=lambda r:r['chi2'])
 uncertainty=local_pair_uncertainty(alt)
 summary={'target':config['target'],'z_abs':config['z'],'lines':config['lines'],'window':config['window'],
   'ncomp':n,'ndata':null['ndata'],'null':null['name'],'alternative':alt['name'],
   'chi2_null':null['chi2'],'chi2_alternative':alt['chi2'],'nominal_ndf_null':null['nominal_ndf'],
   'delta_chi2':null['chi2']-alt['chi2'],'extra_shifts':alt['npar']-null['npar'],
   'both_optimizer_success':bool(null['optimizer_success'] and alt['optimizer_success']),
   'ordinary_gate':ordinary['passes_conditional_adequacy_gate'],'local_pair_uncertainty':uncertainty,
   'null_stationarity':null['stationarity'],'alternative_stationarity':alt['stationarity'],
   'interpretation':'Exploratory, conditional relative-line comparison. No calibrated discovery significance, alpha or cosmic-time-law fit.'}
 (folder/'selected_pair.json').write_text(json.dumps(summary,indent=2)+'\n')
 sig=uncertainty['rank_sensitivity'][1]['conditional_sigma_m_s']
 if sig is not None:
  center=uncertainty['shift_m_s']/1000;values=sorted(set(float(np.clip(center+k*sig/1000,-.999,.999)) for k in [-2,-1,0,1,2]))
  profiles=[]
  for i,x in enumerate(values):
   prof=fit(config,n,True,alt['seed'],alt,name=f'n{n}_profile2382_{i}_os21',oversample=21,max_nfev=300,fixed_shift=x)
   profiles.append({'value_m_s':1000*x,'chi2':prof['chi2'],'delta_from_free':prof['chi2']-alt['chi2'],'optimizer_success':prof['optimizer_success'],'name':prof['name']})
  (folder/'pair_profile.json').write_text(json.dumps({'free_name':alt['name'],'free_chi2':alt['chi2'],'grid':profiles,'interpretation':'Sparse nuisance-reoptimized profile diagnostic; not guaranteed confidence interval; report any lower profile point as an optimization warning.'},indent=2)+'\n')
 return summary

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('configuration');parser.add_argument('--workers',type=int,default=3);args=parser.parse_args()
 print(json.dumps(run_campaign(args.configuration,args.workers),indent=2))
