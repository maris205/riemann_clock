#!/usr/bin/env python3
"""Independent non-optimizing review of combined core masking and empirical GLS.

Reconstruct masks without binary_dilation, construct covariance independently,
replay outputs and selection. Writes only new joint_controls_review artifacts.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,datetime
import numpy as np
from espresso_conventional_null import Model,ROOT
from espresso_saturation_controls import CoreMaskModel
from joint_controls_fit import JointModel
OUT=ROOT/'results/joint_controls_review';OUT.mkdir(parents=True,exist_ok=True)
FIT=ROOT/'results/joint_controls'
checks=[];cases=[];replays=[];comparisons=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,passed,**details):checks.append(dict(name=name,passed=bool(passed),**details))
def same(name,a,b,tol=1e-9):
 a=np.asarray(a);b=np.asarray(b)
 error=float(np.max(np.abs(a-b))) if a.size and a.shape==b.shape else (0. if a.shape==b.shape else float('inf'))
 check(name,a.shape==b.shape and error<=tol,max_abs_difference=error,tolerance=tol)
def covariance(native,kernel):
 # Independent row-by-row construction; do not reuse fit factors or matrix.
 answer=np.zeros((len(native),len(native)))
 for i,pixel in enumerate(native):
  for lag,value in enumerate(kernel):
   answer[i,np.abs(native-pixel)==lag]=value
 return answer
controls=json.loads((ROOT/'results/noise_covariance/controls.json').read_text())
acf=np.array(controls['primary']['acf']);kernel=np.array([x*(len(acf)-i)/len(acf) for i,x in enumerate(acf)])
variance=controls['primary']['variance']
source_paths=[ROOT/'code/joint_controls_fit.py',ROOT/'code/espresso_saturation_controls.py',ROOT/'code/espresso_conventional_null.py',ROOT/'results/espresso_null/null_cross.npz',ROOT/'results/noise_covariance/controls.json',Path(__file__)]
sources={str(p.relative_to(ROOT)):sha(p) for p in source_paths}
prior=json.loads((ROOT/'results/followup_review/validation.json').read_text())
for p in ['code/espresso_saturation_controls.py','code/espresso_conventional_null.py','results/espresso_null/null_cross.npz','results/noise_covariance/controls.json']:
 check('historical_source_unchanged_'+p,sources[p]==prior['source_hashes'][p],sha256=sources[p])
historical_paths=[ROOT/'results/followup_review/validation.json',ROOT/'results/saturation_controls/core01_comparison.json',ROOT/'results/saturation_controls/core03_padded_comparison.json',ROOT/'results/empirical_gls/comparison.json']
historical_before={str(p.relative_to(ROOT)):sha(p) for p in historical_paths}
frozen=np.load(ROOT/'results/espresso_null/null_cross.npz')
full=Model(free_shifts=True)
prior_null=json.loads((ROOT/'results/espresso_null/null_cross.json').read_text())
prior_alt=json.loads((ROOT/'results/espresso_null/alternative_from_published.json').read_text())
p0=np.array(prior_null['parameters']);p1=np.array(prior_alt['parameters'])
fr,fj,fp=full.evaluate(p1)
models={}
for name,threshold,padding in [('core01_gls',.1,0.),('core03_gls',.3,1.)]:
 null=JointModel(threshold=threshold,padding_km_s=padding)
 alt=JointModel(threshold=threshold,padding_km_s=padding,free_shifts=True)
 raw=CoreMaskModel(threshold=threshold,padding_km_s=padding,free_shifts=True)
 models[(name,False)]=null;models[(name,True)]=alt
 check(name+'_nested_labels',null.labels==alt.labels[:null.npar] and alt.npar-null.npar==5)
 check(name+'_expected_dimensions',null.ndata==(2692 if threshold==.1 else 2473) and null.npar==147 and alt.npar==152)
 same(name+'_kernel',alt.kernel,kernel,1e-15)
 rn,jn,_=null.evaluate(p0);ra,ja,_=alt.evaluate(np.r_[p0,np.zeros(5)])
 same(name+'_zero_shift_nesting_residual',rn,ra,0.)
 same(name+'_zero_shift_nesting_jacobian',jn,ja[:,:null.npar],0.)
 rr,jr,rawprofiles=raw.evaluate(p1);rw,jw,profiles=alt.evaluate(p1)
 totalq=0.;gradient=np.zeros(alt.npar);blocks=[];counts=[];gap_details=[]
 for i,(L,N,B,sl,bsl,chol) in enumerate(zip(alt.lines,null.lines,full.lines,alt.slices,full.slices,alt.chols)):
  k=L['key'];remove=np.zeros(len(B['good']),bool);npad=int(np.ceil(padding/full.dv))
  if k in [2382,2600]:
   centers=np.flatnonzero(frozen[f'{k}_model']<threshold)
   if len(centers):remove=np.min(np.abs(np.arange(len(remove))[:,None]-centers[None,:]),axis=1)<=npad
  expected=B['good']&~remove
  check(f'{name}_{k}_independent_frozen_mask',np.array_equal(expected,L['good']))
  check(f'{name}_{k}_H0_H1_same_mask',np.array_equal(N['good'],L['good']))
  check(f'{name}_{k}_no_added_pixels',np.all(~L['good']|B['good']))
  if k not in [2382,2600]:check(f'{name}_{k}_weak_mask_unchanged',np.array_equal(L['good'],B['good']))
  same(f'{name}_{k}_same_native_indices',L['idx'],B['idx'],0.)
  keep=L['good'][B['good']]
  same(f'{name}_{k}_raw_residual_subset',rr[sl],fr[bsl][keep],1e-10)
  same(f'{name}_{k}_raw_jacobian_subset',jr[sl],fj[bsl][keep],1e-10)
  same(f'{name}_{k}_unchanged_profile',profiles[i],fp[i],0.)
  native=L['idx'][L['good']];cov=covariance(native,kernel)
  fullcov=covariance(B['idx'][B['good']],kernel)
  same(f'{name}_{k}_principal_submatrix',cov,fullcov[np.ix_(keep,keep)],0.)
  same(f'{name}_{k}_Cholesky',chol@chol.T,cov,1e-12)
  same(f'{name}_{k}_same_pair_covariance',null.chols[i]@null.chols[i].T,cov,1e-12)
  same(f'{name}_{k}_unit_diagonal',np.diag(cov),np.ones(len(native)),0.)
  eig=np.linalg.eigvalsh(cov)
  check(f'{name}_{k}_positive_definite',eig[0]>0,min_eigenvalue=float(eig[0]))
  solved=np.linalg.solve(cov,rr[sl]);q=float(rr[sl]@solved);totalq+=q;gradient+=jr[sl].T@solved
  same(f'{name}_{k}_dense_quadratic',[rw[sl]@rw[sl]],[q],1e-8)
  for at in np.flatnonzero(np.diff(native)>1):
   lag=int(native[at+1]-native[at]);expected_c=float(kernel[lag]) if lag<len(kernel) else 0.
   same(f'{name}_{k}_gap_at_{int(at)}',[cov[at,at+1]],[expected_c],0.)
   if lag<len(kernel):check(f'{name}_{k}_short_gap_not_reset_{int(at)}',cov[at,at+1]!=0.)
   gap_details.append(dict(line=k,left_retained_row=int(at),native_distance=lag,covariance=expected_c,within_kernel=lag<len(kernel)))
  blocks.append(dict(line=k,npix=len(native),min_eigenvalue=float(eig[0]),max_eigenvalue=float(eig[-1])))
  count=dict(line=k,original=int(B['good'].sum()),removed=int((B['good']&remove).sum()),retained=int(expected.sum()))
  audit=alt.mask_audit[i]
  check(f'{name}_{k}_audit_counts',audit['original_pixels']==count['original'] and audit['removed_pixels']==count['removed'] and audit['retained_pixels']==count['retained'])
  counts.append(count)
 same(name+'_total_dense_quadratic',[rw@rw],[totalq],1e-8)
 same(name+'_dense_gradient',jw.T@rw,gradient,1e-7)
 # All model parameter families, every additional shift, both masked cases.
 labels=['logN_0','logN_22','v_0','v_22','logb_0','logb_22','continuum0_2260','continuum1_2600']+[f'shift_{k}' for k in alt.shift_keys]
 fd=[]
 for label in labels:
  at=alt.labels.index(label);step=2e-4 if label.startswith(('v_','shift_')) else 1e-5
  if label.startswith('continuum'):step=1e-6
  up=p1.copy();down=p1.copy();up[at]+=step;down[at]-=step
  derivative=(alt.evaluate(up)[0]-alt.evaluate(down)[0])/(2*step)
  rel=float(np.linalg.norm(derivative-jw[:,at])/max(np.linalg.norm(jw[:,at]),1e-12))
  check(name+'_finite_difference_'+label,rel<2e-5,relative_l2_error=rel,step=step)
  fd.append(dict(label=label,relative_l2_error=rel,step=step))
 rfirst,jfirst,_=alt.evaluate(p1);rrepeat,jrepeat,_=alt.evaluate(p1.copy())
 check(name+'_cache_same_parameter_values',rfirst is rrepeat and jfirst is jrepeat)
 mutation=p1.copy();mutation[-1]+=.02;alt.evaluate(mutation);mutation[-1]+=.013
 changed=alt.evaluate(mutation);fresh=JointModel(threshold=threshold,padding_km_s=padding,free_shifts=True)
 independently=fresh.evaluate(mutation)
 same(name+'_cache_mutation_residual',changed[0],independently[0],0.)
 same(name+'_cache_mutation_jacobian',changed[1],independently[1],0.)
 roundtrip=alt.evaluate(p1)
 same(name+'_cache_roundtrip_residual',roundtrip[0],rw,0.)
 same(name+'_cache_roundtrip_jacobian',roundtrip[1],jw,0.)
 check(name+'_fit_does_not_move_mask',all(np.array_equal(L['good'],N['good']) for L,N in zip(alt.lines,null.lines)))
 cases.append(dict(case=name,threshold=threshold,requested_padding_km_s=padding,padding_native_pixels=npad,effective_padding_km_s=npad*full.dv,ndata=alt.ndata,mask_counts=counts,covariance_blocks=blocks,gaps=gap_details,finite_differences=fd))
# Replay every complete, final or unfinished optimization endpoint currently on disk.
for path in sorted(FIT.glob('*.json')):
 d=json.loads(path.read_text())
 if 'optimizer_success' not in d:continue
 npz=path.with_suffix('.npz')
 if not npz.exists():continue
 case='core03_gls' if path.stem.startswith('core03_gls') else 'core01_gls'
 m=models[(case,d['free_shifts'])];p=np.array(d['parameters']);r,J,profiles=m.evaluate(p)
 a=np.load(npz);prefix=path.stem
 check(prefix+'_labels',d['labels']==m.labels)
 check(prefix+'_dimensions',d['ndata']==m.ndata and d['npar']==m.npar)
 check(prefix+'_configuration',d['threshold']==(.3 if case=='core03_gls' else .1) and d['padding_km_s']==(1. if case=='core03_gls' else 0.))
 same(prefix+'_chi2',[r@r],[d['chi2']],1e-8)
 same(prefix+'_saved_residual',r,a['residuals'],1e-10)
 same(prefix+'_saved_jacobian',J,a['jacobian'],1e-9)
 same(prefix+'_saved_parameters',p,a['parameters'],0.)
 same(prefix+'_kernel',d['covariance_kernel'],kernel,1e-15)
 same(prefix+'_continuum_variance',[d['continuum_variance']],[variance],0.)
 initial,lo,hi=m.initial();check(prefix+'_bounds',np.all(p>=lo-1e-10) and np.all(p<=hi+1e-10))
 g=J.T@r;distance=np.ones_like(p);distance[g<0]=hi[g<0]-p[g<0];distance[g>0]=p[g>0]-lo[g>0]
 same(prefix+'_optimality',[np.max(np.abs(g*distance))],[d['optimality']],1e-6)
 check(prefix+'_status_matches_success',d['optimizer_success']==(d['status']>0))
 source=json.loads((ROOT/d['start_file']).read_text());lookup=dict(zip(source['labels'],source['parameters']))
 mapped=np.clip([lookup.get(k,v) for k,v in zip(m.labels,initial)],lo+1e-9,hi-1e-9)
 same(prefix+'_mapped_start',mapped,d['initial_parameters'],0.)
 for file,hashvalue in d['source_sha256'].items():check(prefix+'_source_hash_'+file,sha(ROOT/file)==hashvalue)
 for L,prof in zip(m.lines,profiles):
  k=L['key'];same(f'{prefix}_{k}_profile',prof,a[f'{k}_model'],1e-11)
  same(f'{prefix}_{k}_native_indices',L['idx'],a[f'{k}_idx'],0.)
  check(f'{prefix}_{k}_mask',np.array_equal(L['good'],a[f'{k}_good']))
 replays.append(dict(name=prefix,case=case,free_shifts=d['free_shifts'],chi2=d['chi2'],optimizer_success=d['optimizer_success'],status=d['status'],message=d['message'],optimality=d['optimality'],active_bound_count=len(d['active_bounds']),nfev=d['nfev'],json_sha256=sha(path),npz_sha256=sha(npz)))
for case in ['core01_gls','core03_gls']:
 path=FIT/(case+'_comparison.json')
 if not path.exists():comparisons.append(dict(case=case,status='PENDING'));continue
 d=json.loads(path.read_text());a=json.loads((FIT/(d['null']+'.json')).read_text());b=json.loads((FIT/(d['alternative']+'.json')).read_text())
 same(case+'_selected_delta',[a['chi2']-b['chi2']],[d['delta_chi2']],1e-10)
 check(case+'_selected_H1_not_worse_than_H0',b['chi2']<=a['chi2']+1e-8)
 check(case+'_selected_counts',a['ndata']==b['ndata']==d['ndata'] and b['npar']-a['npar']==d['added_parameters']==5)
 check(case+'_selected_flags',d['null_success']==a['optimizer_success'] and d['alternative_success']==b['optimizer_success'])
 trials=[x for x in replays if x['case']==case]
 check(case+'_all_finished_trials_recorded',set(x['name'] for x in trials)==set(x['name'] for x in d['recorded_trials']))
 for recorded in d['recorded_trials']:
  actual=next((x for x in trials if x['name']==recorded['name']),None)
  check(case+'_recorded_trial_'+recorded['name'],actual is not None and all(actual[k]==recorded[k] for k in ['name','free_shifts','chi2','optimizer_success','nfev','optimality']))
 for free,selected in [(False,a),(True,b)]:
  available=[x for x in trials if x['free_shifts']==free]
  same(case+('_H1' if free else '_H0')+'_minimum_ALL_endpoints',[selected['chi2']],[min(x['chi2'] for x in available)],1e-10)
  check(case+('_H1' if free else '_H0')+'_no_omitted_lower_capped',not any(x['chi2']<selected['chi2']-1e-8 for x in available if not x['optimizer_success']))
 # Explicitly evaluate the embedded H1 at the final selected H0 endpoint.
 n=models[(case,False)];alt=models[(case,True)];pn=np.array(a['parameters'])
 nr,nj,_=n.evaluate(pn);ar,aj,_=alt.evaluate(np.r_[pn,np.zeros(5)])
 same(case+'_FINAL_H0_embedded_H1_residual',nr,ar,0.)
 same(case+'_FINAL_H0_embedded_H1_jacobian',nj,aj[:,:n.npar],0.)
 for k in [2260,2344,2374,2382,2586,2600]:
  an=np.load(FIT/(d['null']+'.npz'));bn=np.load(FIT/(d['alternative']+'.npz'))
  check(f'{case}_{k}_FINAL_same_pixels',np.array_equal(an[f'{k}_good'],bn[f'{k}_good']))
 # Recompute both candidate H0 initial scores, then verify chosen initial source.
 scores=[]
 for candidate in d['initial_null_candidates']:
  old=json.loads((ROOT/candidate['file']).read_text());lookup=dict(zip(old['labels'],old['parameters']))
  initial,lo,hi=n.initial();initial=np.clip([lookup.get(k,v) for k,v in zip(n.labels,initial)],lo+1e-9,hi-1e-9)
  rr=n.evaluate(initial)[0];score=float(rr@rr);scores.append(score)
  same(case+'_initial_score_'+candidate['file'],[score],[candidate['chi2']],1e-8)
 initial_fit=json.loads((FIT/(case+'_null.json')).read_text())
 check(case+'_initial_null_minimum_source',initial_fit['start_file']==d['initial_null_candidates'][int(np.argmin(scores))]['file'])
 comparisons.append(dict(status='REPLAYED',summary_sha256=sha(path),selected_both_terminated=a['optimizer_success'] and b['optimizer_success'],**d))
for p,oldhash in historical_before.items():check('historical_file_not_modified_'+p,sha(ROOT/p)==oldhash)
result=dict(verification='PASS' if all(x['passed'] for x in checks) else 'FAIL',timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_checks=len(checks),n_pass=sum(x['passed'] for x in checks),n_fail=sum(not x['passed'] for x in checks),all_comparisons_completed=all(x['status']=='REPLAYED' for x in comparisons),all_selected_optimizers_terminated=all(x.get('selected_both_terminated',False) for x in comparisons),source_sha256=sources,historical_artifact_sha256=historical_before,raw_rho1=float(acf[1]),tapered_kernel=kernel.tolist(),continuum_variance=variance,cases=cases,fit_replays=replays,comparisons=comparisons,checks=checks,scope='Independent mask/covariance/derivative/cache/selection and saved-objective replay only; no optimization rerun or global optimum, calibrated significance, alpha or cosmic-time-law certification.')
(OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['verification','n_checks','n_pass','n_fail','all_comparisons_completed','all_selected_optimizers_terminated']}))
for c in checks:
 if not c['passed']:print('FAIL',c)
print(json.dumps([{k:x.get(k) for k in ['case','status','null','alternative','delta_chi2','selected_both_terminated']} for x in comparisons],indent=2))
