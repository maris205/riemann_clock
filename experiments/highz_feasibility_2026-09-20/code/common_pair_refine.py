#!/usr/bin/env python3
"""Preserved-basin companion: fair H0/H1 and common-pair profile refinement.

The first-pass numerical products remain unchanged. Its final summary exporter
had a duplicate metadata key after all fits were saved. This separate entry
point reads those saved products, checks the improved nuisance basin, and
writes new refined products plus the final uncertainty summary.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import json,time,datetime
import numpy as np
from scipy.optimize import least_squares
import common_pair_fit as base
from archive_profile_expected_noise import ExpectedModel
ROOT=base.ROOT;OUT=base.OUT;SOURCE=base.SOURCE
read=base.read;write=base.write;sha=base.sha

def path_record(path):
 r=read(path);r['fit_file']=str(Path(path).relative_to(ROOT))
 if base.KEY in r.get('labels',[]):r['D_m_s']=1000*r['parameters'][r['labels'].index(base.KEY)]
 return r

def fit_h0(name,model,start,meta,max_nfev=400):
 folder=OUT/model.target;path=folder/(name+'.json')
 if path.exists():
  r=read(path)
  if r['source_sha256']!=meta['source_sha256'] or r['start_sha256']!=sha(start):raise RuntimeError('Existing H0 differs')
  return r
 old=read(start);lookup=dict(zip(old['labels'],old['parameters']));p,lo,hi=model.initial();p=np.clip([lookup.get(k,v) for k,v in zip(model.labels,p)],lo+1e-9,hi-1e-9);initial=p.copy();began=time.time()
 opt=least_squares(model.fun,p,jac=model.jac,bounds=(lo,hi),x_scale='jac',ftol=2e-8,xtol=1e-8,gtol=1e-6,max_nfev=max_nfev)
 r,J,profiles=model.evaluate(opt.x)
 result={**meta,'name':name,'fit_file':str(path.relative_to(ROOT)),'hypothesis':'all_relative_line_shifts_fixed_zero','target':model.target,'z_abs':model.cfg['z'],'ncomp':model.ncomp,'labels':model.labels,'parameters':opt.x.tolist(),'initial_parameters':initial.tolist(),'bounds_lower':lo.tolist(),'bounds_upper':hi.tolist(),'optimized_parameter_indices':list(range(model.npar)),'fixed_shift_km_s':None,'D_m_s':0.,'ndata':model.ndata,'npar_total':model.npar,'npar_optimized':model.npar,'chi2':float(r@r),'optimizer_success':bool(opt.success),'optimizer_status':int(opt.status),'message':opt.message,'nfev':int(opt.nfev),'max_nfev':max_nfev,'optimality':float(opt.optimality),'elapsed_seconds':time.time()-began,'active_bounds':[model.labels[i] for i in range(model.npar) if min(opt.x[i]-lo[i],hi[i]-opt.x[i])<1e-5],'shifts_m_s':{},'start_file':str(Path(start).relative_to(ROOT)),'start_sha256':sha(start)}
 write(path,result);arrays=dict(parameters=opt.x,residual=r,jacobian=J,optimized_parameter_indices=np.arange(model.npar))
 for L,prof in zip(model.lines,profiles):
  for key in ['v','wave','flux','error','good','source_indices']:arrays[f'{L["key"]}_{key}']=L[key]
  arrays[f'{L["key"]}_model']=prof
 np.savez_compressed(path.with_suffix('.npz'),**arrays)
 print('DONE',model.target,name,result['chi2'],opt.success,opt.nfev,flush=True);return result

def run_target(selected):
 target=selected['target'];folder=OUT/target
 model=ExpectedModel(target,selected['ncomp'],True,33);nullmodel=ExpectedModel(target,selected['ncomp'],False,33)
 meta=base.metadata(model,selected);meta['source_sha256'][str(Path(__file__).relative_to(ROOT))]=sha(__file__)
 meta['first_pass_export_note']='First-pass summary construction raised duplicate selection metadata key after numeric products were saved; this companion assembles final summaries from preserved products.'
 oldalt=path_record(SOURCE/(selected['expected_alternative']+'.json'));oldnull=path_record(SOURCE/(selected['expected_null']+'.json'))
 first=[]
 for p in sorted(folder.glob('*.json')):
  if p.name.startswith('refined_'):continue
  r=read(p)
  if 'parameters' in r and 'optimizer_success' in r:first.append(path_record(p))
 if not first:raise RuntimeError('First-pass numeric endpoints required')
 first_hashes={r['fit_file']:sha(ROOT/r['fit_file']) for r in first}
 free=[oldalt]+[r for r in first if r['fixed_shift_km_s'] is None];allh1=[oldalt]+first
 reference=min(allh1,key=lambda r:r['chi2'])
 nulls=[oldnull,fit_h0('refined_h0_from_best_h1',nullmodel,ROOT/reference['fit_file'],meta)]
 h0=min(nulls,key=lambda r:r['chi2']);free.append(base.fit('refined_h1_from_h0',model,ROOT/h0['fit_file'],None,meta,400));allh1.append(free[-1]);reference=min(allh1,key=lambda r:r['chi2'])
 initial_reference=reference;local=base.local_information(model,reference['parameters'])
 _,lo,hi=model.initial();at=model.labels.index(base.KEY)
 points=[];records=[];counter=0
 def selected_existing(theta):return [r for r in allh1 if r.get('fixed_shift_km_s') is not None and abs(r['fixed_shift_km_s']-theta)<1e-12]
 def do_point(theta,phase):
  nonlocal counter
  if any(abs(p['fixed_shift_km_s']-theta)<1e-12 for p in points):return
  name=f'refined_profile_{counter:03d}_{phase}';counter+=1
  warm=min([reference]+points,key=lambda r:abs(r['D_m_s']/1000-theta))
  a=base.fit(name+'_warm',model,ROOT/warm['fit_file'],float(theta),meta)
  b=base.fit(name+'_cross_h0',model,ROOT/h0['fit_file'],float(theta),meta)
  records.extend([a,b]);allh1.extend([a,b]);points.append(min(selected_existing(theta),key=lambda r:r['chi2']))
  write(folder/'refined_profile_progress.json',dict(target=target,points=[{k:r[k] for k in ['fit_file','D_m_s','chi2','optimizer_success']} for r in sorted(points,key=lambda r:r['D_m_s'])]))
 if target=='J232128-105122':
  # The first-pass two-start profile did not reveal a distinct basin here.
  coordinates=sorted(set(r['fixed_shift_km_s'] for r in first if r['fixed_shift_km_s'] is not None))
  points=[min(selected_existing(theta),key=lambda r:r['chi2']) for theta in coordinates]
 else:
  center=reference['D_m_s']/1000;sig=local['conditional_sigma_m_s']/1000
  grid=np.unique(np.r_[np.clip(center+np.arange(-3,3.01,.5)*sig,lo[at],hi[at]),center,0.])
  write(folder/'refined_design.json',dict(target=target,grid_center_km_s=center,local_sigma_km_s=sig,grid_km_s=grid.tolist(),source_reference=reference['fit_file'],other_floating_line_shifts=[k for k in model.shiftkeys if k!=2382],fixed_line_shift=2382,shift_bounds_km_s=[float(lo[at]),float(hi[at])],max_nfev_per_fit=300,starts_per_point=2,**meta))
  for theta in sorted(grid,key=lambda x:abs(x-center)):do_point(float(theta),'grid')
  reference=min(allh1,key=lambda r:r['chi2'])
  bestfree=min(free,key=lambda r:r['chi2'])
  if reference['chi2']<bestfree['chi2']-1e-5:
   free.append(base.fit('refined_h1_from_grid_minimum',model,ROOT/reference['fit_file'],None,meta,400));allh1.append(free[-1]);reference=min(allh1,key=lambda r:r['chi2'])
  for level in [1.,3.84]:
   for iteration in range(2):
    contours=base.crossings(points,reference,level)
    for side in ['lower','upper']:
     b=contours[side]
     if b:do_point(float(np.mean(b['bracket_m_s'])/1000),f'refine{level:g}_{side}_{iteration}')
    reference=min(allh1,key=lambda r:r['chi2'])
  # One bounded opposite traversal checks initialization hysteresis at the
  # already specified grid; no new masks/components/shift bounds are selected.
  for i,p in enumerate(sorted(points,key=lambda r:r['D_m_s'],reverse=True)):
   theta=p['fixed_shift_km_s'];start=min(allh1,key=lambda r:r['chi2'])
   r=base.fit(f'refined_reverse_{i:03d}',model,ROOT/start['fit_file'],theta,meta)
   records.append(r);allh1.append(r)
  points=[min(selected_existing(theta),key=lambda r:r['chi2']) for theta in sorted(set(p['fixed_shift_km_s'] for p in points))]
 reference=min(allh1,key=lambda r:r['chi2']);bestfree=min(free,key=lambda r:r['chi2'])
 if reference['chi2']<bestfree['chi2']-1e-5:
  free.append(base.fit('refined_h1_from_final_profile',model,ROOT/reference['fit_file'],None,meta,400));allh1.append(free[-1]);reference=min(allh1,key=lambda r:r['chi2'])
 # Give the fully conventional comparison the same improved nuisance basin.
 if reference['chi2']<initial_reference['chi2']-1e-5:
  nulls.append(fit_h0('refined_h0_from_final_h1',nullmodel,ROOT/reference['fit_file'],meta))
  h0=min(nulls,key=lambda r:r['chi2'])
  free.append(base.fit('refined_h1_from_final_h0',model,ROOT/h0['fit_file'],None,meta,400));allh1.append(free[-1]);reference=min(allh1,key=lambda r:r['chi2'])
 h0=min(nulls,key=lambda r:r['chi2']);bestfree=min(free,key=lambda r:r['chi2'])
 local=base.local_information(model,reference['parameters']);write(folder/'refined_local_information.json',dict(**local,reference_file=reference['fit_file'],**meta))
 intervals=[base.crossings(points,reference,level) for level in [1.,3.84]]
 selected_points=[dict(D_m_s=r['D_m_s'],chi2=r['chi2'],delta_chi2=r['chi2']-reference['chi2'],fit_file=r['fit_file'],optimizer_success=r['optimizer_success'],active_bounds=r['active_bounds'],nfev=r['nfev'],optimality=r['optimality']) for r in sorted(points,key=lambda r:r['D_m_s'])]
 profile_gap=bestfree['chi2']-reference['chi2'];resolved=profile_gap<=1e-4
 first_vs_refined=[]
 for p in points:
  old=[r for r in first if r['fixed_shift_km_s'] is not None and abs(r['fixed_shift_km_s']-p['fixed_shift_km_s'])<1e-12]
  if old:first_vs_refined.append(dict(D_m_s=p['D_m_s'],first_chi2=min(r['chi2'] for r in old),refined_chi2=p['chi2']))
 summary={**meta,'target':target,'z_abs':model.cfg['z'],'D_m_s':reference['D_m_s'],'conditional_sigma_m_s':local['conditional_sigma_m_s'],'point_estimate_reference_file':reference['fit_file'],'point_estimate_reference_sha256':sha(ROOT/reference['fit_file']),'reference_chi2':reference['chi2'],'reference_optimizer_success':reference['optimizer_success'],'reference_is_unrestricted':reference in free,'unrestricted_reference_gap_chi2':profile_gap,'profile_baseline_resolved_to_tolerance':resolved,'baseline_tolerance_chi2':1e-4,'profile_grid_D_m_s':[r['D_m_s'] for r in selected_points],'profile_delta_chi2':[r['delta_chi2'] for r in selected_points],'profile_points':selected_points,'profile_intervals':intervals,'active_bounds':reference['active_bounds'],'profile_all_selected_terminated':all(r['optimizer_success'] for r in points),'profile_contours_closed':all(x['closed'] for x in intervals),'profile_grid_reaches_shift_bound':bool(any(abs(r['D_m_s']/1000-lo[at])<1e-9 or abs(r['D_m_s']/1000-hi[at])<1e-9 for r in points)),'ndata':model.ndata,'ncomp':model.ncomp,'npar':model.npar,'sourcefit_paths':[oldnull['fit_file'],oldalt['fit_file']],'sourcefit_hashes':{r['fit_file']:sha(ROOT/r['fit_file']) for r in [oldnull,oldalt]},'original_D_m_s':oldalt['D_m_s'],'original_h1_chi2':oldalt['chi2'],'original_h0_chi2':oldnull['chi2'],'full_comparison':dict(null_file=h0['fit_file'],alternative_file=reference['fit_file'],chi2_null=h0['chi2'],chi2_alternative=reference['chi2'],delta_chi2=h0['chi2']-reference['chi2'],added_parameters=model.npar-nullmodel.npar,null_optimizer_success=h0['optimizer_success'],alternative_optimizer_success=reference['optimizer_success'],null_active_bounds=h0['active_bounds'],null_npar=nullmodel.npar,alternative_npar=model.npar,interpretation='All original relative line offsets fixed zero for H0. D=0 profile leaves other offsets free and is a different model.'),'all_profile_fit_files':sorted(set(r['fit_file'] for r in first+records if r.get('fixed_shift_km_s') is not None)),'full_fit_files':[r['fit_file'] for r in free],'null_fit_files':[r['fit_file'] for r in nulls],'first_pass_sha256':first_hashes,'first_vs_refined_at_common_grid':first_vs_refined,'endpoint_selection':'Lowest recorded objective including unfinished endpoints; profile grid retains minimum across all starts at each final coordinate.','interpretation':'Conditional fixed-architecture common-pair likelihood contour in preferred row2 diagonal errors. No calibration/covariance/model uncertainty marginalization, no coverage-calibrated interval or physical alpha/time measurement.'}
 write(folder/'refined_summary.json',summary)
 assert all(sha(ROOT/f)==h for f,h in first_hashes.items()),'First-pass products changed'
 return summary

def main():
 selected=read(SOURCE/'summary.json')['comparisons']
 with ProcessPoolExecutor(2) as pool:rows=list(pool.map(run_target,selected))
 summary=dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),common_pair='D=v(FeII2382)-v(FeII2374), m/s',comparisons=rows,source_stage='refined common-pair nuisance profile; historical row2 results preserved',interpretation='Descriptive conditional profile intervals; not a cosmic-time-law inference.')
 write(OUT/'refined_summary.json',summary)
 print(json.dumps([{k:r[k] for k in ['target','D_m_s','conditional_sigma_m_s','reference_chi2','profile_contours_closed','profile_all_selected_terminated','profile_baseline_resolved_to_tolerance','full_comparison']} for r in rows],indent=2),flush=True)
if __name__=='__main__':main()
