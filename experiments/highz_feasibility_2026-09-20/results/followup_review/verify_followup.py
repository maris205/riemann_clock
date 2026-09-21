#!/usr/bin/env python3
"""Independent, non-optimizing verification of saturation and empirical GLS fits."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import sys,json,hashlib,datetime
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'code'))
from espresso_conventional_null import Model
from espresso_saturation_controls import CoreMaskModel
from espresso_empirical_gls import GLSModel
OUT=ROOT/'results/followup_review'
checks=[];metrics={};replays=[];comparisons=[]
def check(name,ok,**details):
 checks.append(dict(name=name,passed=bool(ok),**details))
def same(name,a,b,tol=1e-9):
 a=np.asarray(a);b=np.asarray(b)
 err=float(np.max(np.abs(a-b))) if a.size else 0.
 check(name,a.shape==b.shape and err<=tol,max_abs_difference=err,tolerance=tol)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
core_hash=sha(ROOT/'code/espresso_conventional_null.py')
check('unchanged_original_physical_core',core_hash=='424dc9e32cb9bb0a3eece533c64d44ae2582f1edc4b8fa4b14889c6e4db142de',sha256=core_hash)
source_hashes={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'code/espresso_saturation_controls.py',ROOT/'code/espresso_empirical_gls.py',ROOT/'code/polish_espresso_followup.py',ROOT/'code/espresso_conventional_null.py',ROOT/'results/espresso_null/null_cross.npz',ROOT/'results/espresso_null/null_cross.json',ROOT/'results/noise_covariance/controls.json',Path(__file__)]}
frozen=np.load(ROOT/'results/espresso_null/null_cross.npz')
base=Model();baseline=json.loads((ROOT/'results/espresso_null/null_cross.json').read_text())
p0=np.array(baseline['parameters']);br,bj,bprofiles=base.evaluate(p0)
for L,prof in zip(base.lines,bprofiles):
 same(f'frozen_null_{L["key"]}_physical_replay',prof,frozen[f'{L["key"]}_model'],1e-11)
 check(f'frozen_null_{L["key"]}_mask_replay',np.array_equal(L['good'],frozen[f'{L["key"]}_good']))
same('frozen_null_source_objective',[br@br],[baseline['chi2']],1e-8)
mask_counts=[]
for case,threshold,padding in [('core01',.1,0.),('core03_padded',.3,1.)]:
 null=CoreMaskModel(threshold=threshold,padding_km_s=padding)
 alt=CoreMaskModel(threshold=threshold,padding_km_s=padding,free_shifts=True)
 check(case+'_parameter_nesting',null.labels==alt.labels[:null.npar] and alt.npar-null.npar==5)
 r,j,_=null.evaluate(p0)
 altp=np.r_[p0,np.zeros(5)]
 altr,altj,_=alt.evaluate(altp)
 same(case+'_zero_shift_nesting_residual',r,altr,0.)
 same(case+'_zero_shift_nesting_nuisance_jacobian',j,altj[:,:null.npar],0.)
 expected_r=[];expected_j=[];computed_counts=[]
 for i,(L,A,B,sl) in enumerate(zip(null.lines,alt.lines,base.lines,base.slices)):
  k=L['key'];remove=np.zeros(len(L['good']),bool);pad=int(np.ceil(padding/base.dv))
  if k in [2382,2600]:
   centers=np.flatnonzero(frozen[f'{k}_model']<threshold)
   # Independent direct finite-distance construction (no binary_dilation).
   if len(centers):remove=np.min(np.abs(np.arange(len(remove))[:,None]-centers[None,:]),axis=1)<=pad
  expected=B['good']&~remove
  check(f'{case}_{k}_independent_mask',np.array_equal(L['good'],expected))
  check(f'{case}_{k}_same_H0_H1_pixels',np.array_equal(L['good'],A['good']))
  check(f'{case}_{k}_same_native_wavelength',np.array_equal(L['wave'],B['wave']))
  check(f'{case}_{k}_no_added_pixels',np.all(~L['good']|B['good']))
  kept=expected[B['good']]
  expected_r.append(br[sl][kept]);expected_j.append(bj[sl][kept])
  count=dict(line=k,original=int(B['good'].sum()),removed=int((B['good']&remove).sum()),retained=int(expected.sum()))
  computed_counts.append(count)
  audit=null.mask_audit[i]
  check(f'{case}_{k}_audit_counts',audit['removed_pixels']==count['removed'] and audit['retained_pixels']==count['retained'])
  if k not in [2382,2600]:check(f'{case}_{k}_weak_line_untouched',np.array_equal(L['good'],B['good']))
 same(case+'_residual_is_fixed_subset',r,np.concatenate(expected_r),1e-10)
 same(case+'_jacobian_is_fixed_subset',j,np.vstack(expected_j),1e-10)
 masks_before=[L['good'].copy() for L in null.lines]
 shifted=p0.copy();shifted[0]+=.01;null.evaluate(shifted)
 check(case+'_mask_does_not_move_with_fit',all(np.array_equal(x,L['good']) for x,L in zip(masks_before,null.lines)))
 mask_counts.append(dict(case=case,threshold=threshold,requested_padding_km_s=padding,dilation_pixels=pad,actual_padding_km_s=pad*base.dv,ndata=null.ndata,per_line=computed_counts))
# Build the covariance independently without using GLSModel's stored factors.
controls=json.loads((ROOT/'results/noise_covariance/controls.json').read_text())
acf=np.asarray(controls['primary']['acf']);kernel=np.array([v*(len(acf)-i)/len(acf) for i,v in enumerate(acf)])
gls=GLSModel(free_shifts=True)
altbaseline=json.loads((ROOT/'results/espresso_null/alternative_from_published.json').read_text())
p=np.array(altbaseline['parameters'])
raw=Model(free_shifts=True)
r0,j0,profiles0=raw.evaluate(p);rw,jw,profiles=gls.evaluate(p)
same('Bartlett_kernel_reconstructed',gls.kernel,kernel,1e-15)
quad=0.;gradient=np.zeros(gls.npar);glsdetails=[]
for L,sl,Lchol,profile,rawprofile in zip(gls.lines,gls.slices,gls.cholesky,profiles,profiles0):
 k=L['key'];native=L['idx'][L['good']];n=len(native)
 cov=np.zeros((n,n))
 for i in range(n):
  separation=np.abs(native[i]-native)
  use=separation<len(kernel);cov[i,use]=kernel[separation[use]]
 same(f'gls_{k}_cholesky_reconstructs_covariance',Lchol@Lchol.T,cov,1e-12)
 eig=np.linalg.eigvalsh(cov)
 check(f'gls_{k}_positive_definite',eig[0]>0,min_eigenvalue=float(eig[0]),max_eigenvalue=float(eig[-1]))
 same(f'gls_{k}_same_forward_spectrum',profile,rawprofile,0.)
 check(f'gls_{k}_unit_diagonal_without_variance_transfer',np.array_equal(np.diag(cov),np.ones(n)))
 solved=np.linalg.solve(cov,r0[sl]);q=float(r0[sl]@solved);quad+=q
 gradient+=j0[sl].T@solved
 same(f'gls_{k}_quadratic_vs_dense_solve',[rw[sl]@rw[sl]],[q],1e-8)
 # Sparse masks must use original index distance, not compressed row index.
 gaps=np.flatnonzero(np.diff(native)>1)
 for at in gaps:
  lag=int(native[at+1]-native[at]);expected=float(kernel[lag]) if lag<len(kernel) else 0.
  same(f'gls_{k}_masked_gap_{int(at)}',[cov[at,at+1]],[expected],0.)
 glsdetails.append(dict(line=k,npix=n,min_eigenvalue=float(eig[0]),max_eigenvalue=float(eig[-1]),masked_gaps=int(len(gaps))))
same('GLS_total_quadratic_independent_dense_solve',[rw@rw],[quad],1e-8)
same('GLS_gradient_independent_dense_solve',jw.T@rw,gradient,1e-7)
glsnull=GLSModel();gnr,gnj,_=glsnull.evaluate(p0)
gnar,gnaj,_=gls.evaluate(np.r_[p0,np.zeros(5)])
same('GLS_zero_shift_nesting_residual',gnr,gnar,0.)
same('GLS_zero_shift_nesting_nuisance_jacobian',gnj,gnaj[:,:glsnull.npar],0.)
# Finite differences at a saved physical alternative, covering every parameter family
# and all five relative-shift derivatives; these test the changed whitening path.
fdlabels=['logN_0','logN_22','v_0','v_22','logb_0','logb_22','continuum0_2260','continuum1_2600']+[f'shift_{k}' for k in gls.shift_keys]
fd=[]
for label in fdlabels:
 i=gls.labels.index(label);step=2e-4 if label.startswith(('v_','shift_')) else 1e-5
 if label.startswith('continuum'):step=1e-6
 up=p.copy();dn=p.copy();up[i]+=step;dn[i]-=step
 rp=gls.evaluate(up)[0];rm=gls.evaluate(dn)[0];numeric=(rp-rm)/(2*step)
 denom=max(np.linalg.norm(jw[:,i]),1e-12);rel=float(np.linalg.norm(numeric-jw[:,i])/denom)
 check('GLS_finite_difference_'+label,rel<2e-5,relative_l2_error=rel,step=step)
 fd.append(dict(label=label,relative_l2_error=rel,step=step))
# Cache returns both consistent results and recomputes after in-place mutation.
ra,ja,pa=gls.evaluate(p);rb,jb,pb=gls.evaluate(p.copy())
check('GLS_cache_repeated_values',ra is rb and ja is jb)
q=p.copy();q[-1]+=.013;rq,jq,_=gls.evaluate(q);q[-1]+=.017
rqq,jqq,_=gls.evaluate(q);fresh=GLSModel(free_shifts=True);rf,jf,_=fresh.evaluate(q)
same('GLS_cache_inplace_parameter_change_residual',rqq,rf,0.)
same('GLS_cache_inplace_parameter_change_jacobian',jqq,jf,0.)
same('GLS_cache_roundtrip_residual',gls.evaluate(p)[0],rw,0.)
same('GLS_cache_roundtrip_jacobian',gls.evaluate(p)[1],jw,0.)
metrics.update(mask_cases=mask_counts,covariance=dict(raw_rho1=float(acf[1]),tapered_rho1=float(kernel[1]),kernel=kernel.tolist(),continuum_variance=controls['primary']['variance'],blocks=glsdetails),finite_differences=fd)
# Saved results are replayed but no optimization is re-run. Pick up every finished file.
for folder in ['saturation_controls','empirical_gls']:
 for path in sorted((ROOT/'results'/folder).glob('*.json')):
  d=json.loads(path.read_text())
  if 'optimizer_success' not in d:continue
  arrpath=path.with_suffix('.npz')
  if not arrpath.exists():continue
  kwargs=dict(free_shifts=d['free_shifts'])
  m=CoreMaskModel(threshold=d['threshold'],padding_km_s=d['padding_km_s'],**kwargs) if folder=='saturation_controls' else GLSModel(**kwargs)
  pp=np.array(d['parameters']);rr,jj,pr=m.evaluate(pp);a=np.load(arrpath);name=folder+'/'+path.stem
  check(name+'_labels',d['labels']==m.labels)
  check(name+'_dimensions',d['ndata']==m.ndata and d['npar']==m.npar and len(rr)==m.ndata)
  same(name+'_chi2_replay',[rr@rr],[d['chi2']],1e-8)
  same(name+'_saved_residuals',rr,a['residuals'],1e-10)
  same(name+'_saved_jacobian',jj,a['jacobian'],1e-9)
  same(name+'_saved_parameters',pp,a['parameters'],0.)
  _,lo,hi=m.initial();check(name+'_parameters_in_bounds',np.all(pp>=lo-1e-10) and np.all(pp<=hi+1e-10))
  gradient=jj.T@rr
  distance=np.ones_like(pp)
  distance[gradient<0]=hi[gradient<0]-pp[gradient<0]
  distance[gradient>0]=pp[gradient>0]-lo[gradient>0]
  replayed_optimality=float(np.max(np.abs(gradient*distance)))
  same(name+'_bounded_gradient_optimality',[replayed_optimality],[d['optimality']],1e-6)
  if folder=='saturation_controls':
   expected_threshold=.3 if path.stem.startswith('core03') else .1
   expected_padding=1. if path.stem.startswith('core03') else 0.
   check(name+'_scenario_configuration',d['threshold']==expected_threshold and d['padding_km_s']==expected_padding)
  start=ROOT/d['start_file'];old=json.loads(start.read_text());lookup=dict(zip(old['labels'],old['parameters']));ip=m.initial()[0]
  ip=np.array([lookup.get(k,v) for k,v in zip(m.labels,ip)]);ip=np.clip(ip,lo+1e-9,hi-1e-9)
  same(name+'_start_mapping_replay',ip,d['initial_parameters'],0.)
  for L,prof in zip(m.lines,pr):
   k=L['key'];same(name+f'_{k}_profile',prof,a[f'{k}_model'],1e-11)
   check(name+f'_{k}_saved_mask',np.array_equal(L['good'],a[f'{k}_good']))
  if folder=='saturation_controls':check(name+'_frozen_mask_hash',d['mask_source_sha256']==sha(ROOT/'results/espresso_null/null_cross.npz'))
  else:
   check(name+'_covariance_source_hash',d['covariance_source_sha256']==sha(ROOT/'results/noise_covariance/controls.json'))
   same(name+'_separate_variance_rescaling',[d['chi2_if_continuum_variance_transferred']],[d['chi2']/d['continuum_variance']],1e-8)
  replays.append(dict(name=name,optimizer_success=d['optimizer_success'],message=d['message'],chi2=d['chi2'],nfev=d['nfev'],optimality=d['optimality'],active_bound_count=len(d['active_bounds']),sha256=sha(path),npz_sha256=sha(arrpath)))
for folder,filename in [('saturation_controls','core01_comparison.json'),('saturation_controls','core03_padded_comparison.json'),('empirical_gls','comparison.json')]:
 path=ROOT/'results'/folder/filename
 if not path.exists():comparisons.append(dict(file=str(path.relative_to(ROOT)),status='PENDING'));continue
 d=json.loads(path.read_text());a=json.loads((path.parent/(d['null']+'.json')).read_text());b=json.loads((path.parent/(d['alternative']+'.json')).read_text())
 name=folder+'/'+filename
 same(name+'_delta_replay',[a['chi2']-b['chi2']],[d['delta_chi2']],1e-10)
 check(name+'_nested_dimensions',a['npar']+5==b['npar'] and a['ndata']==b['ndata']==d['ndata'])
 check(name+'_both_optimizer_success',a['optimizer_success'] and b['optimizer_success'])
 # Audit all completed starts of this scenario, including unsuccessful ones.
 attempted=[]
 prefix=d.get('case','tapered')
 for candidate in path.parent.glob(prefix+'*.json'):
  c=json.loads(candidate.read_text())
  if 'optimizer_success' not in c:continue
  attempted.append(c)
 for free,selected in [(False,a),(True,b)]:
  candidates=[c for c in attempted if c['free_shifts']==free]
  successful=[c for c in candidates if c['optimizer_success']]
  if successful:
   same(name+('_H1' if free else '_H0')+'_best_tested_successful',[selected['chi2']],[min(c['chi2'] for c in successful)],1e-10)
  unfinished_lower=[dict(name=c['name'],chi2=c['chi2']) for c in candidates if not c['optimizer_success'] and c['chi2']<selected['chi2']-1e-6]
  check(name+('_H1' if free else '_H0')+'_no_lower_unterminated_fit',not unfinished_lower,lower_unterminated_fits=unfinished_lower)
 for k in [2260,2344,2374,2382,2586,2600]:
  an=np.load(path.parent/(d['null']+'.npz'));bn=np.load(path.parent/(d['alternative']+'.npz'))
  check(name+f'_{k}_identical_pair_mask',np.array_equal(an[f'{k}_good'],bn[f'{k}_good']))
 if 'delta_chi2_if_continuum_variance_transferred' in d:same(name+'_delta_variance_rescaling',[d['delta_chi2_if_continuum_variance_transferred']],[d['delta_chi2']/d['continuum_variance']],1e-10)
 comparisons.append(dict(file=str(path.relative_to(ROOT)),status='REPLAYED',sha256=sha(path),**d))
result=dict(verification='PASS' if all(c['passed'] for c in checks) else 'FAIL',all_comparisons_completed=all(x['status']=='REPLAYED' for x in comparisons),timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='Independent fixed-mask reconstruction, covariance algebra/finite derivatives/cache checks, source hashes and saved fit replay; no optimization rerun/global optimum certification.',n_checks=len(checks),n_pass=sum(c['passed'] for c in checks),n_fail=sum(not c['passed'] for c in checks),source_hashes=source_hashes,metrics=metrics,fit_replays=replays,comparisons=comparisons,checks=checks)
(OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['verification','all_comparisons_completed','n_checks','n_pass','n_fail']}))
for c in checks:
 if not c['passed']:print('FAIL',c)
print(json.dumps(comparisons,indent=2))
