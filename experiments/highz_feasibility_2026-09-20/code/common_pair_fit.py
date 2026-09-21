#!/usr/bin/env python3
"""Conditional 2382-minus-2374 uncertainty for two frozen archive pilots.

Preferred expected-fluctuation row2; all original H1 nuisance parameters and
other line shifts float while profiling only shift_2382. New products only.
No covariance/calibration/time-law inference or coverage-calibrated interval.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import argparse,datetime,hashlib,json,time
import numpy as np
from scipy.optimize import least_squares
import archive_profile_pilot as pilot
from archive_profile_expected_noise import ExpectedModel,OUT as SOURCE
ROOT=pilot.ROOT;OUT=ROOT/'results/common_pair';KEY='shift_2382';OS=33
MAX_NFEV=300

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def local_information(model,parameters):
 p=np.asarray(parameters);r,J,_=model.evaluate(p);i=model.labels.index(KEY)
 A=np.delete(J,i,axis=1);scl=np.linalg.norm(A,axis=0)
 if np.any(scl<=0):raise ValueError('Zero nuisance Jacobian column')
 U,s,V=np.linalg.svd(A/scl,full_matrices=False);g=J[:,i];rows=[]
 for cutoff in [1e-8,1e-10,1e-12]:
  keep=s>s[0]*cutoff;projected=g-U[:,keep]@(U[:,keep].T@g);information=float(projected@projected)
  rows.append(dict(relative_svd_cutoff=cutoff,nuisance_rank=int(keep.sum()),nuisance_columns=A.shape[1],projected_information_per_km_s_squared=information,conditional_sigma_m_s=1000/np.sqrt(information),projected_column=projected.tolist()))
 _,lo,hi=model.initial();active=[model.labels[j] for j in range(model.npar) if min(p[j]-lo[j],hi[j]-p[j])<1e-5]
 return dict(target=model.target,D_m_s=float(1000*p[i]),conditional_sigma_m_s=rows[1]['conditional_sigma_m_s'],cutoff_sensitivity=rows,nuisance_normalized_singular_values=s.tolist(),active_bounds=active,chi2=float(r@r),ndata=model.ndata,npar=model.npar,interpretation='Unconstrained local nuisance tangent projection in row2 diagonal weights, no reduced-chi-square rescaling; boundary/model/calibration uncertainty not included.')

def metadata(model,selected):
 files=[Path(__file__),Path(pilot.__file__),Path(__file__).with_name('archive_profile_expected_noise.py'),pilot.ATOMIC,pilot.DATA/(model.target+'.npz'),pilot.DATA/(model.target+'_metadata.json'),SOURCE/'summary.json',SOURCE/(selected['expected_null']+'.json'),SOURCE/(selected['expected_alternative']+'.json')]
 files += [pilot.DATA/f'{model.target}_z{model.cfg["z"]:.3f}_FeII{k}_window.npz' for k in model.keys]
 return dict(source_sha256={str(p.relative_to(ROOT)):sha(p) for p in files},configuration=model.cfg,lines=model.keys,anchor_line=2374,target_line=2382,oversample=OS,error_kind='FITS primary row2 expected fluctuation; diagonal errors without reduced-chi-square rescaling',selection='Frozen row1-selected architecture/windows/pixels; row2 selected H1 nuisance freedom; no model reselection',atomic_response='D=v2382-v2374; sign-equivalent log energy/frequency-ratio perturbation Y=-D/c is kinematic only, not a theory response.')

def fit(name,model,start_file,fixed_shift_km_s,meta,max_nfev=MAX_NFEV):
 folder=OUT/model.target;folder.mkdir(parents=True,exist_ok=True);path=folder/(name+'.json')
 if path.exists():
  previous=read(path)
  if previous['source_sha256']!=meta['source_sha256'] or previous['fixed_shift_km_s']!=fixed_shift_km_s or previous['start_sha256']!=sha(start_file):raise RuntimeError('Existing product differs from requested run: '+str(path))
  if not path.with_suffix('.npz').exists():raise RuntimeError('Incomplete product: '+str(path))
  return previous
 source=read(start_file);lookup=dict(zip(source['labels'],source['parameters']));p,lo,hi=model.initial()
 p=np.clip([lookup.get(k,v) for k,v in zip(model.labels,p)],lo+1e-9,hi-1e-9);at=model.labels.index(KEY)
 free=np.ones(model.npar,bool)
 if fixed_shift_km_s is not None:
  if not lo[at]<=fixed_shift_km_s<=hi[at]:raise ValueError('Profile coordinate outside original bounds')
  p[at]=fixed_shift_km_s;free[at]=False
 initial=p.copy();calls=0;began=time.time()
 def unpack(x):
  answer=initial.copy();answer[free]=x;return answer
 def fun(x):
  nonlocal calls
  pp=unpack(x);r=model.evaluate(pp)[0];calls+=1
  if calls%100==0:
   print(model.target,name,calls,float(r@r),flush=True)
   write(folder/(name+'_checkpoint.json'),dict(labels=model.labels,parameters=pp.tolist(),chi2=float(r@r),calls=calls))
  return r
 def jac(x):return model.evaluate(unpack(x))[1][:,free]
 opt=least_squares(fun,p[free],jac=jac,bounds=(lo[free],hi[free]),x_scale='jac',ftol=2e-8,xtol=1e-8,gtol=1e-6,max_nfev=max_nfev)
 fullp=unpack(opt.x);r,J,profiles=model.evaluate(fullp);chi=float(r@r)
 result=dict(name=name,fit_file=str(path.relative_to(ROOT)),target=model.target,z_abs=model.cfg['z'],ncomp=model.ncomp,labels=model.labels,parameters=fullp.tolist(),initial_parameters=initial.tolist(),bounds_lower=lo.tolist(),bounds_upper=hi.tolist(),optimized_parameter_indices=np.flatnonzero(free).tolist(),fixed_shift_km_s=fixed_shift_km_s,D_m_s=float(1000*fullp[at]),ndata=model.ndata,npar_total=model.npar,npar_optimized=int(free.sum()),chi2=chi,optimizer_success=bool(opt.success),optimizer_status=int(opt.status),message=opt.message,nfev=int(opt.nfev),max_nfev=max_nfev,optimality=float(opt.optimality),elapsed_seconds=time.time()-began,active_bounds=[model.labels[j] for j in range(model.npar) if min(fullp[j]-lo[j],hi[j]-fullp[j])<1e-5],shifts_m_s={str(k):1000*float(fullp[model.labels.index(f'shift_{k}')]) for k in model.shiftkeys},start_file=str(Path(start_file).relative_to(ROOT)),start_sha256=sha(start_file),**meta)
 write(path,result)
 arrays=dict(parameters=fullp,residual=r,jacobian=J,optimized_parameter_indices=np.flatnonzero(free))
 for L,profile in zip(model.lines,profiles):
  for k in ['v','wave','flux','error','good','source_indices']:arrays[f'{L["key"]}_{k}']=L[k]
  arrays[f'{L["key"]}_model']=profile
 np.savez_compressed(path.with_suffix('.npz'),**arrays)
 print('DONE',model.target,name,'D',result['D_m_s'],'chi',chi,'ok',opt.success,'n',opt.nfev,flush=True)
 return result

def crossings(points,reference,level):
 # Return brackets/linear interpolation, never claim calibrated coverage.
 points=sorted(points,key=lambda x:x['D_m_s']);center=reference['D_m_s'];refq=reference['chi2'];sides={}
 for side,sign in [('lower',-1),('upper',1)]:
  inside=dict(D_m_s=center,chi2=refq,optimizer_success=reference['optimizer_success'],fit_file=reference['fit_file'])
  candidates=sorted([p for p in points if sign*(p['D_m_s']-center)>1e-7],key=lambda p:abs(p['D_m_s']-center))
  result=None
  for outside in candidates:
   dy=outside['chi2']-refq
   if dy>=level:
    da=inside['chi2']-refq;db=dy;x=inside['D_m_s']+(outside['D_m_s']-inside['D_m_s'])*(level-da)/(db-da)
    result=dict(crossing_m_s=float(x),bracket_m_s=sorted([inside['D_m_s'],outside['D_m_s']]),bracket_delta_chi2=[da,db],both_endpoints_terminated=inside['optimizer_success'] and outside['optimizer_success'],fit_files=[inside['fit_file'],outside['fit_file']]);break
   inside=outside
  sides[side]=result
 return dict(delta_chi2_level=level,lower=sides['lower'],upper=sides['upper'],closed=sides['lower'] is not None and sides['upper'] is not None,interpretation='Linear interpolation within actual profiled grid brackets; descriptive likelihood-objective contour, not coverage-calibrated interval.')

def run_target(selected):
 target=selected['target'];folder=OUT/target;folder.mkdir(parents=True,exist_ok=True)
 model=ExpectedModel(target,selected['ncomp'],True,OS);meta=metadata(model,selected)
 oldnull=SOURCE/(selected['expected_null']+'.json');oldalt=SOURCE/(selected['expected_alternative']+'.json')
 original=read(oldalt);original['fit_file']=str(oldalt.relative_to(ROOT));original['D_m_s']=1000*original['parameters'][model.labels.index(KEY)]
 original_local=local_information(model,original['parameters'])
 write(folder/'original_local_information.json',dict(**original_local,source_fit=str(oldalt.relative_to(ROOT)),source_fit_sha256=sha(oldalt),**meta))
 full=[original,fit('full_from_h1',model,oldalt,None,meta),fit('full_from_h0',model,oldnull,None,meta)]
 reference=min(full,key=lambda x:x['chi2'])
 local=local_information(model,reference['parameters']);center=reference['D_m_s']/1000;sig=local['conditional_sigma_m_s']/1000
 _,lo,hi=model.initial();at=model.labels.index(KEY)
 grid=np.unique(np.r_[np.clip(center+np.arange(-3,3.01,.5)*sig,lo[at],hi[at]),center,0.])
 design=dict(target=target,z_abs=model.cfg['z'],ncomp=model.ncomp,grid_center_km_s=center,grid_scale_local_sigma_km_s=sig,initial_grid_km_s=grid.tolist(),profile_coordinate=KEY,other_relative_shifts_float=model.shiftkeys,one_coordinate_fixed=True,max_nfev_per_fit=MAX_NFEV,starts_per_point=2,adaptive_rounds_per_contour=2,contours=[1.,3.84],shift_bounds_km_s=[float(lo[at]),float(hi[at])],initial_baseline=reference['fit_file'],**meta)
 write(folder/'design.json',design)
 records=[];points=[];counter=0
 def do_point(theta,phase):
  nonlocal counter
  if any(abs(p['fixed_shift_km_s']-theta)<1e-12 for p in points):return
  name=f'profile_{counter:03d}_{phase}';counter+=1
  # Warm neighboring point and independent null-derived nuisance branch.
  candidates=[reference]+points
  warm=min(candidates,key=lambda x:abs(x['D_m_s']/1000-theta))
  a=fit(name+'_warm',model,ROOT/warm['fit_file'],float(theta),meta)
  b=fit(name+'_cross_null',model,oldnull,float(theta),meta)
  records.extend([a,b]);chosen=min([a,b],key=lambda x:x['chi2'])
  points.append(chosen)
  write(folder/'profile_progress.json',dict(target=target,points=[{k:p[k] for k in ['fit_file','fixed_shift_km_s','D_m_s','chi2','optimizer_success','nfev','optimality']} for p in sorted(points,key=lambda x:x['D_m_s'])],recorded_fits=[r['fit_file'] for r in records]))
 for theta in sorted(grid,key=lambda x:abs(x-center)):do_point(float(theta),'grid')
 # A conditional endpoint below the unrestricted endpoint triggers one recorded
 # full-dimensional restart, with no alterations to the old published pilots.
 lowest=min(points,key=lambda x:x['chi2'])
 if lowest['chi2']<reference['chi2']-1e-5:
  full.append(fit('full_from_profile_minimum',model,ROOT/lowest['fit_file'],None,meta,400))
 reference=min(full+points,key=lambda x:x['chi2'])
 for level in [1.,3.84]:
  for iteration in range(2):
   contours=crossings(points,reference,level)
   for side in ['lower','upper']:
    bracket=contours[side]
    if bracket is not None:do_point(float(np.mean(bracket['bracket_m_s'])/1000),f'refine{level:g}_{side}_{iteration}')
   reference=min(full+points,key=lambda x:x['chi2'])
 lowest=min(points,key=lambda x:x['chi2']);bestfull=min(full,key=lambda x:x['chi2'])
 if lowest['chi2']<bestfull['chi2']-1e-5:
  full.append(fit('full_from_final_profile_minimum',model,ROOT/lowest['fit_file'],None,meta,400))
 reference=min(full+points,key=lambda x:x['chi2']);bestfull=min(full,key=lambda x:x['chi2'])
 local=local_information(model,reference['parameters']);write(folder/'local_information.json',dict(**local,reference_file=reference['fit_file'],**meta))
 points=sorted(points,key=lambda x:x['D_m_s']);intervals=[crossings(points,reference,level) for level in [1.,3.84]]
 # Record every endpoint, not just successful or favorable trials.
 selected_fixed=[dict(D_m_s=p['D_m_s'],chi2=p['chi2'],delta_chi2=p['chi2']-reference['chi2'],fit_file=p['fit_file'],optimizer_success=p['optimizer_success'],active_bounds=p['active_bounds'],nfev=p['nfev'],optimality=p['optimality']) for p in points]
 summary=dict(target=target,z_abs=model.cfg['z'],D_m_s=reference['D_m_s'],conditional_sigma_m_s=local['conditional_sigma_m_s'],point_estimate_reference_file=reference['fit_file'],point_estimate_reference_sha256=sha(ROOT/reference['fit_file']),reference_chi2=reference['chi2'],reference_optimizer_success=reference['optimizer_success'],reference_is_unrestricted=reference in full,unrestricted_reference_gap_chi2=bestfull['chi2']-reference['chi2'],profile_grid_D_m_s=[p['D_m_s'] for p in points],profile_delta_chi2=[p['chi2']-reference['chi2'] for p in points],profile_points=selected_fixed,profile_intervals=intervals,active_bounds=reference['active_bounds'],profile_all_selected_terminated=all(p['optimizer_success'] for p in points),profile_contours_closed=all(x['closed'] for x in intervals),profile_grid_reaches_shift_bound=bool(any(abs(p['D_m_s']/1000-lo[at])<1e-9 or abs(p['D_m_s']/1000-hi[at])<1e-9 for p in points)),ndata=model.ndata,ncomp=model.ncomp,npar=model.npar,sourcefit_paths=[str(oldnull.relative_to(ROOT)),str(oldalt.relative_to(ROOT))],sourcefit_hashes={str(p.relative_to(ROOT)):sha(p) for p in [oldnull,oldalt]},full_fit_files=[r['fit_file'] for r in full],all_profile_fit_files=[r['fit_file'] for r in records],selection='Minimum objective among all recorded endpoints including unfinished points; profile minima use both starts at each fixed coordinate.',inference_scope='Conditional fixed-architecture common-pair uncertainty using row2 diagonal errors. Reference contour levels 1 and 3.84 are not calibrated coverage probabilities. Gas architecture, pixel covariance, residual calibration and cross-target systematics are not marginalized. No physical alpha/time inference.',**meta)
 write(folder/'summary.json',summary)
 return summary

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--target',choices=list(pilot.CONFIG));args=parser.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 selected=read(SOURCE/'summary.json')['comparisons'];selected=[x for x in selected if not args.target or x['target']==args.target]
 if len(selected)==1:rows=[run_target(selected[0])]
 else:
  with ProcessPoolExecutor(2) as pool:rows=list(pool.map(run_target,selected))
 summary=dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),common_pair='D = v(Fe II 2382)-v(Fe II 2374), m/s',comparisons=rows,interpretation='Actual common observable and conditional nuisance-profile uncertainty; calibration/model covariance unknown, no time-law measurement.')
 write(OUT/('summary.json' if not args.target else args.target+'_summary.json'),summary)
 print(json.dumps([{k:r[k] for k in ['target','D_m_s','conditional_sigma_m_s','profile_contours_closed','profile_all_selected_terminated','reference_is_unrestricted']} for r in rows],indent=2),flush=True)
if __name__=='__main__':main()
