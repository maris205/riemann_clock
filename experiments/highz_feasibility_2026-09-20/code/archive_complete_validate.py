#!/usr/bin/env python3
"""Independent numerical/provenance audit of the bounded archive campaign.

Run again after additional targets finish. Prior audit artifacts are untouched.
"""
import argparse
import hashlib
import io
import json
import math
import re
import sys
import tarfile
from pathlib import Path
import numpy as np
from scipy.special import voigt_profile

BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'.deps'))
from astropy.io import fits
import archive_complete_model as core

OUT=core.OUT/'validation'
CHECKS=[]
WARNINGS=[]
FD=[]
QUADRATURE=[]
RAW={}
ACTUAL_STARTS=[]
COUNT_SELECTION=[]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check(name,ok,detail=None):CHECKS.append(dict(name=name,passed=bool(ok),detail=detail))
def close(name,a,b,rtol=3e-10,atol=3e-10):
 a,b=np.asarray(a),np.asarray(b)
 ok=a.shape==b.shape and np.allclose(a,b,rtol=rtol,atol=atol,equal_nan=True)
 detail=float(np.nanmax(np.abs(a.astype(float)-b.astype(float)))) if a.shape==b.shape and a.size else None
 check(name,ok,detail)

def raw_data(target):
 if target in RAW:return RAW[target]
 meta=json.loads((core.DATA/f'{target}_metadata.json').read_text());path=BASE/meta['archive_file']
 check(target+' original archive SHA256',sha(path)==meta['archive_sha256'])
 with tarfile.open(path) as archive:blob=archive.extractfile(f'{target}/{target}.fits').read()
 check(target+' original FITS SHA256',hashlib.sha256(blob).hexdigest()==meta['fits_sha256'])
 with fits.open(io.BytesIO(blob)) as hdus:
  data=hdus[0].data.copy();header=hdus[0].header.copy()
 wave=np.power(10.,header['CRVAL1']+(np.arange(data.shape[1])+1-header['CRPIX1'])*header['CD1_1'])
 valid=(data[4]==1)&np.isfinite(data[0])&np.isfinite(data[1])&(data[1]>0)
 RAW[target]=(data,header,wave,valid,meta)
 return RAW[target]

def atom_rows(key):
 rows=[]
 for line in core.old.ATOMIC.read_text().splitlines():
  fields=line.split()
  if len(fields)>=6 and fields[0]=='FeII' and key<=float(fields[1])<key+1:rows.append([float(x) for x in fields[1:6]])
 return np.asarray(rows)

def model_class(target):
 # The extended J2331 region needs its explicit neighbouring-transition model;
 # do not silently validate it with a single-transition opacity model.
 if target=='J233156-090802':
  import archive_complete_j2331_model as neighbor
  return neighbor.NeighborModel
 return core.Model

def independent_profile(model,p,line):
 n=model.ncomp;grid=line['grid'];key=line['key'];shift=p[model.labels.index(f'shift_{key}')] if key in model.shiftkeys else 0.
 rows=atom_rows(key)
 if model.target=='J233156-090802' and key==1608:rows=np.vstack([rows,atom_rows(1611)])
 # Gaussian cgs e, electron mass and c; coefficient converted from cm and cm/s
 # to Angstrom and km/s. No model optical-depth coefficient is reused.
 e=1.602176634e-19*2.99792458e9;me=9.1093837139e-28;c=2.99792458e10
 coefficient=math.sqrt(math.pi)*e*e/(me*c)*1e-13
 opacity=np.zeros(len(grid))
 for logn,v,lb in zip(p[:n],p[n:2*n],p[2*n:3*n]):
  b=math.exp(lb)
  for wavelength,oscillator,gamma,_mass,_q in rows:
   u=-core.C*np.expm1((v+core.C*math.log(wavelength/line['ref'])+shift-grid)/core.C)/b
   damping=gamma*wavelength*1e-13/(4*math.pi*b)
   opacity+=coefficient*10.**logn*oscillator*wavelength/b*math.sqrt(math.pi)*voigt_profile(u,1/math.sqrt(2),damping)
 i=3*n+4*model.keys.index(key);c0,c1,zero,lf=p[i:i+4]
 sigma=math.exp(lf)/2.354820045/(line['dv']/model.oversample);radius=int(6*sigma+.5)
 x=np.arange(-radius,radius+1);kernel=np.exp(-.5*(x/sigma)**2);kernel/=kernel.sum()
 smooth=np.convolve(np.pad(np.exp(-opacity),(radius,radius),mode='edge'),kernel,mode='valid')
 flux=smooth.reshape(-1,model.oversample).mean(axis=1)[line['pad']:-line['pad']]
 return (1+c0+c1*line['x'])*(zero+(1-zero)*flux)

def audit_data_model(config,cls):
 target=config['target'];data,header,wave,valid,meta=raw_data(target)
 m=cls(config,3,True,9);null=cls(config,3,False,9);p,lo,hi=m.initial(0)
 # Deterministic interior perturbations exercise every derivative family.
 p=np.clip(p,lo+.002,hi-.002)
 for j,label in enumerate(m.labels):
  if label.startswith('shift_'):p[j]=.2
 indices=[]
 for line in m.lines:
  key=line['key'];name=f'{target} FeII{key}';rows=atom_rows(key);ref=np.average(rows[:,0],weights=rows[:,1]);v=core.C*np.log(wave/(ref*(1+config['z'])))
  ids=np.flatnonzero((v>=config['window'][0])&(v<=config['window'][1]))
  close(name+' original native source indices',line['source_indices'],ids,rtol=0,atol=0)
  close(name+' wavelengths reconstructed from original FITS',line['wave'],wave[ids],atol=3e-11)
  close(name+' original valid status mask',line['good'],valid[ids],rtol=0,atol=0)
  close(name+' original expected error row2',line['error'],data[2,ids],rtol=0,atol=0)
  close(name+' original statistical error row1',line['statistical_error'],data[1,ids],rtol=0,atol=0)
  close(name+' original flux',line['flux'],data[0,ids],rtol=0,atol=0)
  close(name+' actual native log dispersion',line['dv'],core.C*math.log(10)*header['CD1_1'])
  check(name+' finite positive expected error on fixed mask',np.all(np.isfinite(line['error'][line['good']])&(line['error'][line['good']]>0)))
  check(name+' contiguous native wavelength window',np.all(np.diff(ids)==1))
  near=min(meta['tables']['1']['rows'],key=lambda x:abs(x['Wavelength']-np.median(wave[ids])))
  close(name+' nominal LSF from coadd metadata',line['nominal_fwhm'],core.C/near['NomResolPower'])
  indices.extend(ids[line['good']].tolist())
  expected_atoms=np.vstack([rows,atom_rows(1611)]) if target=='J233156-090802' and key==1608 else rows
  close(name+' atomic opacity terms without repeated abundance',line['atom'],expected_atoms)
 check(target+' globally unique source indices in likelihood',len(indices)==len(set(indices)),dict(retained=len(indices),unique=len(set(indices))))
 if target=='J233156-090802':check(target+' no separate1611 likelihood region',1611 not in m.keys and 1608 in m.keys)
 n=m.ncomp
 close(target+' logN lower bounds',lo[:n],np.full(n,9.));close(target+' logN upper bounds',hi[:n],np.full(n,17.))
 close(target+' component velocity lower',lo[n:2*n],np.full(n,config['component_bounds'][0]));close(target+' component velocity upper',hi[n:2*n],np.full(n,config['component_bounds'][1]))
 close(target+' logb lower',lo[2*n:3*n],np.full(n,np.log(.5)));close(target+' logb upper',hi[2*n:3*n],np.full(n,np.log(30.)))
 for i,line in enumerate(m.lines):
  a=3*n+4*i;close(target+' line nuisance lower '+str(line['key']),lo[a:a+4],[-.10,-.03,-.02,np.log(.6*line['nominal_fwhm'])]);close(target+' line nuisance upper '+str(line['key']),hi[a:a+4],[.10,.03,.02,np.log(1.4*line['nominal_fwhm'])])
 close(target+' shift lower',lo[-len(m.shiftkeys):],np.full(len(m.shiftkeys),-1.));close(target+' shift upper',hi[-len(m.shiftkeys):],np.full(len(m.shiftkeys),1.))
 r,J,profiles=m.evaluate(p);copy=p.copy();p[0]+=.01;check(target+' cache copies parameter key',not np.array_equal(m.cache[0],p));p[:]=copy
 check(target+' cache same parameter hit',m.evaluate(p)[0] is r)
 for j,label in enumerate(m.labels):
  h=.001 if label.startswith(('v_','shift_')) else (1e-5 if label.startswith(('cont','zero')) else 1e-4)
  plus,minus=p.copy(),p.copy();plus[j]+=h;minus[j]-=h
  fd=(m.fun(plus)-m.fun(minus))/(2*h);err=np.linalg.norm(fd-J[:,j])/max(np.linalg.norm(J[:,j]),1e-4)
  FD.append(dict(target=target,label=label,relative_l2_error=float(err)))
  check(target+' finite difference '+label,err<8e-5,float(err))
 for line,profile in zip(m.lines,profiles):close(target+' independent cgs Voigt convolution '+str(line['key']),profile,independent_profile(m,p,line),rtol=2e-8,atol=2e-8)
 pnull=p[:null.npar];pzero=p.copy();pzero[null.npar:]=0
 close(target+' exact nested null profile',m.fun(pzero),null.fun(pnull),rtol=0,atol=0)
 for a,b in zip(m.lines,null.lines):close(target+' identical H0 H1 retained indices '+str(a['key']),a['source_indices'][a['good']],b['source_indices'][b['good']],rtol=0,atol=0)

def audit_fit(path,cls):
 obj=json.loads(path.read_text());target=obj['target'];m=cls(obj['configuration'],obj['ncomp'],obj['free_shifts'],obj['oversample']);p=np.array(obj['parameters']);r,J,profiles=m.evaluate(p);prefix=target+'/'+obj['name'];saved=np.load(path.with_suffix('.npz'))
 check(prefix+' labels',obj['labels']==m.labels)
 active=np.ones(m.npar,bool)
 if obj['fixed_2382_shift_km_s'] is not None:
  j=m.labels.index('shift_2382');active[j]=False;close(prefix+' fixed profile shift',p[j],obj['fixed_2382_shift_km_s'])
  actual=np.asarray(obj['initial_parameters']).copy();actual[j]=obj['fixed_2382_shift_km_s']
  ACTUAL_STARTS.append(dict(target=target,fit_name=obj['name'],source_json=str(path.relative_to(BASE)),source_json_sha256=sha(path),labels=obj['labels'],recorded_initial_parameters=obj['initial_parameters'],fixed_coordinate='shift_2382',fixed_value_km_s=obj['fixed_2382_shift_km_s'],actual_initial_parameters=actual.tolist(),interpretation='Supplemental reconstruction only: overwrite the fixed coordinate after the recorded seed snapshot. Active optimized coordinates and saved fit results are unmodified.'))
  WARNINGS.append(dict(kind='profile_initial_vector_provenance',fit=prefix,detail='Recorded initial_parameters is the free seed before fixed2382 insertion. Actual full initial vector is reconstructed by replacing shift_2382 with fixed_2382_shift_km_s; active optimizer coordinates are unchanged.'))
 check(prefix+' dimensions and fixed parameter accounting',obj['ndata']==m.ndata and obj['npar']==int(active.sum()) and obj['nominal_ndf']==m.ndata-int(active.sum()))
 close(prefix+' fitted parameter mask',saved['fitted_parameter_mask'],active,rtol=0,atol=0)
 close(prefix+' residual replay',r,saved['residual']);close(prefix+' Jacobian replay',J,saved['jacobian']);close(prefix+' parameter replay',p,saved['parameters'],rtol=0,atol=0)
 close(prefix+' objective replay',r@r,obj['chi2']);close(prefix+' reduced objective',obj['chi2_per_ndf'],(r@r)/obj['nominal_ndf'])
 k=int(active.sum());close(prefix+' AICc arithmetic',obj['aicc'],r@r+2*k+2*k*(k+1)/(m.ndata-k-1))
 close(prefix+' model lower bounds',obj['bounds_lower'],m.initial(obj['seed'])[1]);close(prefix+' model upper bounds',obj['bounds_upper'],m.initial(obj['seed'])[2])
 check(prefix+' feasible parameters',np.all(p>=obj['bounds_lower']) and np.all(p<=obj['bounds_upper']))
 check(prefix+' stop flag fidelity',obj['optimizer_success']==(obj['optimizer_status']>0))
 source=obj['source_hashes'];check(prefix+' coadd hash',source['coadd']==sha(core.DATA/f'{target}.npz'));check(prefix+' metadata hash',source['metadata']==sha(core.DATA/f'{target}_metadata.json'));check(prefix+' atomic hash',source['atomic']==sha(core.old.ATOMIC));check(prefix+' base hash',source['base_model']==sha(core.old.__file__))
 check(prefix+' adapter hash',source['adapter']==sha(core.__file__))
 if target=='J233156-090802':
  import archive_complete_j2331_model as neighbor
  check(prefix+' neighbor adapter hash in frozen configuration',obj['configuration']['neighbor_opacity_adapter_sha256']==sha(neighbor.__file__))
 indices=[]
 for line,profile,sl,stated in zip(m.lines,profiles,m.slices,obj['per_line']):
  name=prefix+' '+str(line['key']);close(name+' saved model',profile,saved[f'{line["key"]}_model']);close(name+' line objective',r[sl]@r[sl],stated['chi2']);close(name+' line objective per pixel',r[sl]@r[sl]/line['good'].sum(),stated['chi2_per_pixel'])
  for key in ['v','wave','flux','error','statistical_error','good','source_indices']:close(name+' saved '+key,line[key],saved[f'{line["key"]}_{key}'],rtol=0,atol=0)
  indices.extend(line['source_indices'][line['good']].tolist())
 check(prefix+' global unique likelihood indices',len(indices)==len(set(indices)))
 for key,shift in obj['shifts_m_s'].items():close(prefix+' shift units '+key,shift,1000*p[m.labels.index('shift_'+key)])
 return obj

def audit_count_selection(target,stage,chosen,options):
 good=[x for x in options if x['optimizer_success']]
 expected=min(good or options,key=lambda r:r['chi2'])
 check(target+' '+stage+' good-or-options selector n'+str(chosen['ncomp']),chosen['name']==expected['name'])
 incomplete=[x for x in options if not x['optimizer_success']]
 lower=min(incomplete,key=lambda r:r['chi2']) if incomplete else None
 same=lower is None or lower['oversample']==chosen['oversample']
 check(target+' '+stage+' count comparison same quadrature n'+str(chosen['ncomp']),same)
 COUNT_SELECTION.append(dict(target=target,stage=stage,ncomp=chosen['ncomp'],selected_name=chosen['name'],selected_chi2=chosen['chi2'],selected_optimizer_success=chosen['optimizer_success'],oversample=chosen['oversample'],successful_candidates=len(good),incomplete_candidates=len(incomplete),actual_rule='minimum chi2 among successful stops if any, otherwise minimum chi2 among incomplete endpoints',lowest_incomplete_name=lower['name'] if lower else None,lowest_incomplete_chi2=lower['chi2'] if lower else None,selected_minus_lowest_incomplete_chi2=chosen['chi2']-lower['chi2'] if lower else None,lower_incomplete_available=bool(lower and lower['chi2']<chosen['chi2']),candidate_names=[r['name'] for r in options]))

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--targets',nargs='*');args=parser.parse_args()
 folders=[p for p in core.OUT.iterdir() if p.is_dir() and p.name.startswith('J') and (not args.targets or p.name in args.targets)]
 counts={};selected=[];skipped=[]
 for folder in sorted(folders):
  paths=[]
  for p in folder.rglob('*.json'):
   if p.with_suffix('.npz').exists():
    o=json.loads(p.read_text())
    if {'configuration','parameters','source_hashes','ncomp'}<=set(o):paths.append(p)
  if not paths:skipped.append(folder.name);continue
  try:cls=model_class(folder.name)
  except ModuleNotFoundError:skipped.append(folder.name);continue
  config=json.loads(paths[0].read_text())['configuration'];audit_data_model(config,cls)
  check(folder.name+' unique saved fit record names',len({p.stem for p in paths})==len(paths))
  records={p.stem:audit_fit(p,cls) for p in sorted(paths)};record_paths={p.stem:p for p in paths};counts[folder.name]=len(records)
  ordinary_path=folder/'ordinary_model_summary.json'
  if not ordinary_path.exists():continue
  ordinary=json.loads(ordinary_path.read_text());refined=records[ordinary['selected_null']]
  gate=bool(refined['optimizer_success'] and refined['chi2_per_ndf']<=1.5 and max(x['chi2_per_pixel'] for x in refined['per_line'])<=1.8)
  check(folder.name+' ordinary adequacy gate arithmetic',ordinary['passes_conditional_adequacy_gate']==gate)
  bycount=ordinary['by_count'];eligible=[r for r in bycount if r['optimizer_success']];chosen=min(eligible or bycount,key=lambda r:r['aicc'])
  check(folder.name+' null-only AICc selection',chosen['ncomp']==ordinary['ncomp']==refined['ncomp'])
  for row in bycount:
   fit=records[row['name']]
   for k in ['ncomp','chi2','nominal_ndf','aicc','optimizer_success']:close(folder.name+' bycount '+row['name']+' '+k,row[k],fit[k])
   pattern=rf'n{row["ncomp"]}_null_(s[0-9]+|continued)_os9'
   choices=[x for x in records.values() if re.fullmatch(pattern,x['name'])]
   if choices:audit_count_selection(folder.name,'original_count_grid',fit,choices)
  original_gate=gate;extension_path=folder/'ordinary_model_extension_summary.json'
  if extension_path.exists():
   extension=json.loads(extension_path.read_text());extended=records[extension['selected_null']]
   egate=bool(extended['optimizer_success'] and extended['chi2_per_ndf']<=1.5 and max(x['chi2_per_pixel'] for x in extended['per_line'])<=1.8)
   check(folder.name+' bounded extension gate arithmetic',extension['passes_conditional_adequacy_gate']==egate)
   choices=[records[x['name']] for x in bycount]
   for item in extension['extension_count_records']:
    fit=records[item['name']];choices.append(fit)
    for key in ['ncomp','chi2','aicc','optimizer_success']:close(folder.name+' bounded extension record '+item['name']+' '+key,item[key],fit[key])
   for number in sorted({r['ncomp'] for r in extension['extension_count_records']}):
    opts=[records[x['name']] for x in extension['extension_count_records'] if x['ncomp']==number];good=[x for x in opts if x['optimizer_success']];chosen_count=min(good or opts,key=lambda x:x['chi2']);audit_count_selection(folder.name,'bounded_extension_count_grid',chosen_count,opts)
   choice=min([r for r in choices if r['optimizer_success']],key=lambda r:r['aicc'])
   check(folder.name+' bounded extension null-only successful AICc selection',choice['ncomp']==extended['ncomp']==extension['ncomp'])
   protocol=json.loads((folder/extension['extension_protocol']).read_text())
   check(folder.name+' extension frozen original configuration hash',protocol['unchanged_configuration_sha256']==sha(folder/'selection_frozen.json'))
   check(folder.name+' extension preserved original summary hash',protocol['original_null_summary_sha256']==sha(ordinary_path))
   ordinary_for_pair=extension;gate=egate
  else:ordinary_for_pair=ordinary
  pair_path=folder/'selected_pair.json';selected_names=[ordinary['selected_null']];diagnostic_names=[];diagnostic_summary=None
  if pair_path.exists():
   pair=json.loads(pair_path.read_text());a,b=records[pair['null']],records[pair['alternative']]
   close(folder.name+' selected delta objective',pair['delta_chi2'],a['chi2']-b['chi2'])
   check(folder.name+' selected pair dimension',a['ndata']==b['ndata']==pair['ndata'] and pair['extra_shifts']==b['npar']-a['npar'])
   check(folder.name+' successful stop summary',pair['both_optimizer_success']==(a['optimizer_success'] and b['optimizer_success']))
   check(folder.name+' pair authorized only after applicable gate',pair['ordinary_gate']==gate==True)
   check(folder.name+' pair retains applicable selected component count',a['ncomp']==b['ncomp']==ordinary_for_pair['ncomp'])
   if 'recovery_protocol' in pair:
    recovery=json.loads((folder/pair['recovery_protocol']).read_text());previous=json.loads((folder/pair['supersedes']).read_text())
    check(folder.name+' recovery preserves previous pair snapshot',recovery['previous_pair']==previous)
    endpoint=records[recovery['source_low_endpoint']]
    check(folder.name+' recovery source endpoint hash',recovery['source_sha256']==sha(record_paths[endpoint['name']]))
    recovered=records['n30_null_recovered_os21']
    expected_start=np.clip(np.asarray(endpoint['parameters']),np.asarray(recovered['bounds_lower'])+1e-8,np.asarray(recovered['bounds_upper'])-1e-8)
    close(folder.name+' recovery actual initial vector from lower endpoint',recovered['initial_parameters'],expected_start,rtol=0,atol=0)
    check(folder.name+' recovery same count/lines/window',recovered['ncomp']==endpoint['ncomp'] and recovered['configuration']==endpoint['configuration'])
   selected_names=[pair['null'],pair['alternative']]
   incomplete=[x for x in records.values() if not x['optimizer_success'] and not x['free_shifts'] and x['ncomp']==a['ncomp'] and x['configuration']['lines']==a['configuration']['lines'] and x['configuration']['window']==a['configuration']['window']]
   if incomplete:
    endpoint=min(incomplete,key=lambda x:x['chi2']);em=cls(endpoint['configuration'],endpoint['ncomp'],False,21);er=em.fun(np.asarray(endpoint['parameters']));objective=float(er@er)
    if objective<min(a['chi2'],b['chi2']):WARNINGS.append(dict(kind='lower_incomplete_conventional_endpoint',target=folder.name,name=endpoint['name'],source_oversample=endpoint['oversample'],comparison_oversample=21,endpoint_chi2_at21=objective,selected_null_chi2=a['chi2'],selected_alternative_chi2=b['chi2'],below_selected_null=a['chi2']-objective,below_selected_alternative=b['chi2']-objective,interpretation='Same selected architecture and same quadrature: an unfinished conventional endpoint fits better than the reported successful pair. The pair delta is conditional on selected optimization basins, not a global alternative preference. No new optimization was performed.'))
  elif not gate:check(folder.name+' no automatic primary H1 after failed gate',not any(r['free_shifts'] and r['configuration']['lines']==config['lines'] and record_paths[r['name']].parent==folder for r in records.values()))
  diagnostic_path=folder/'failure_pair/summary.json'
  if diagnostic_path.exists():
   raw_ds=json.loads(diagnostic_path.read_text());ds=dict(raw_ds)
   if isinstance(ds['null'],dict):
    dn0,da0=ds['null'],ds['alternative']
    ds=dict(null=dn0['name'],alternative=da0['name'],chi2_null=dn0['chi2'],chi2_alternative=da0['chi2'],delta_chi2=raw_ds['delta_chi2'],primary_remains_gated_out=bool(raw_ds['primary_adequacy_gate']==False and raw_ds['diagnostic_only'] and raw_ds['primary_precision_descriptor_D'].startswith('not estimated')),ndata=raw_ds['native_pixels'],both_optimizer_success=raw_ds['both_selected_optimizer_success'],extra_shift_parameters=raw_ds['extra_region_shifts'],alternative_passes_same_residual_gate=raw_ds['alternative_meets_same_numeric_quality_thresholds'])
   dn,da=records[ds['null']],records[ds['alternative']]
   close(folder.name+' explanatory failure-pair null objective',ds['chi2_null'],dn['chi2']);close(folder.name+' explanatory failure-pair alternative objective',ds['chi2_alternative'],da['chi2']);close(folder.name+' explanatory failure-pair delta',ds['delta_chi2'],dn['chi2']-da['chi2'])
   check(folder.name+' failure-pair retains failed-primary classification',ds['primary_remains_gated_out']==True and original_gate==False)
   check(folder.name+' failure-pair same full-primary architecture',dn['ncomp']==da['ncomp']==refined['ncomp'] and dn['configuration']['lines']==da['configuration']['lines']==refined['configuration']['lines'] and dn['configuration']['window']==da['configuration']['window']==refined['configuration']['window'])
   check(folder.name+' failure-pair fixed data size',dn['ndata']==da['ndata']==refined['ndata']==ds['ndata'])
   check(folder.name+' failure-pair stop flag fidelity',ds['both_optimizer_success']==(dn['optimizer_success'] and da['optimizer_success']))
   check(folder.name+' failure-pair dimension',ds['extra_shift_parameters']==da['npar']-dn['npar'])
   for fit,key in [(dn,'null_passes_original_conditional_gate'),(da,'alternative_passes_same_residual_gate')]:
    fg=bool(fit['optimizer_success'] and fit['chi2_per_ndf']<=1.5 and max(x['chi2_per_pixel'] for x in fit['per_line'])<=1.8)
    if key in ds:check(folder.name+' failure-pair '+key,ds[key]==fg)
   protocol_path=folder/'failure_pair/protocol_amendment.json'
   if not protocol_path.exists():protocol_path=folder/'failure_pair/protocol.json'
   amendment=json.loads(protocol_path.read_text())
   engine_hash=amendment.get('engine_sha256',amendment.get('shared_engine_sha256'))
   check(folder.name+' failure-pair declared engine hash',engine_hash==sha(core.__file__))
   if 'original_primary_sha256' in amendment:check(folder.name+' failure-pair protected primary source',amendment['original_primary_sha256']==sha(folder/(ds['original_primary_null']+'.json')))
   if 'initial_primary_null_sha256' in amendment:check(folder.name+' failure-pair protected primary source',amendment['initial_primary_null_sha256']==sha(BASE/amendment['initial_primary_null_file']))
   if 'protocol_sha256' in raw_ds:check(folder.name+' failure-pair protocol hash',raw_ds['protocol_sha256']==sha(protocol_path))
   diagnostic_names=[ds['null'],ds['alternative']];diagnostic_summary=dict(null=ds['null'],alternative=ds['alternative'],delta_chi2=ds['delta_chi2'],primary_remains_gated_out=ds['primary_remains_gated_out'])
  metrics=[]
  for name in dict.fromkeys(selected_names+diagnostic_names):
   obj=records[name];p=np.array(obj['parameters']);models=[cls(obj['configuration'],obj['ncomp'],obj['free_shifts'],s) for s in [21,41]];res=[m.fun(p) for m in models];q=dict(target=folder.name,name=name,chi2_21=float(res[0]@res[0]),chi2_41=float(res[1]@res[1]),max_profile_change_in_error_units=float(np.max(np.abs(res[1]-res[0]))));q['delta_chi2_41_minus21']=q['chi2_41']-q['chi2_21'];QUADRATURE.append(q);metrics.append(q)
   check(folder.name+' selected quadrature '+name,abs(q['delta_chi2_41_minus21'])<.1 and q['max_profile_change_in_error_units']<.02,q)
  selected.append(dict(target=folder.name,original_passes_conditional_adequacy_gate=original_gate,passes_conditional_adequacy_gate=gate,selected_null=ordinary_for_pair['selected_null'],selected_names=selected_names,explanatory_failure_diagnostic=diagnostic_summary,quadrature=metrics))
 OUT.mkdir(parents=True,exist_ok=True)
 result=dict(scope='Independent raw-FITS/native-pixel uniqueness, expected errors, adapter bounds/Jacobian, cgs forward profiles, fit replay, selection/gate arithmetic and quadrature. Successful numerical stop is not stationarity or physical-model certification.',passed=all(c['passed'] for c in CHECKS),n_checks=len(CHECKS),n_saved_fits=sum(counts.values()),fits_per_target=counts,targets_without_final_fits=skipped,validator_sha256=sha(__file__),adapter_sha256=sha(core.__file__),max_finite_difference_relative_l2_error=max((x['relative_l2_error'] for x in FD),default=None),selected=selected,quadrature=QUADRATURE,warnings=WARNINGS,failures=[c for c in CHECKS if not c['passed']],checks=CHECKS)
 (OUT/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n')
 (OUT/'profile_actual_initial_parameters.json').write_text(json.dumps(dict(scope='Supplemental reconstruction of actual full starting vectors for fixed2382 profiles; original engine and results are preserved.',profiles=ACTUAL_STARTS),indent=2)+'\n')
 (OUT/'per_count_selection_audit.json').write_text(json.dumps(dict(scope='Actual good-or-options count selection, retaining every lower unfinished endpoint. These gaps describe stage-specific optimizer bookkeeping, not global model ranking; later recovery fits are separate products.',counts=COUNT_SELECTION,lower_incomplete_gaps=[x for x in COUNT_SELECTION if x['lower_incomplete_available']]),indent=2)+'\n')
 lines=['# Independent bounded archive campaign validation','',f'Passed: **{result["passed"]}**. {len(CHECKS)} checks; {result["n_saved_fits"]} stored fits; targets: {", ".join(counts)}.','',f'No final fits yet checked for: {", ".join(skipped) or "none"}.','', 'The audit reconstructs original native wavelength indices and all likelihood arrays from archived FITS, explicitly rejects repeated global native pixels across likelihood regions, verifies frozen H0/H1 pixels and bounds, finite-differences every parameter family, and computes representative flux profiles independently with a cgs optical-depth coefficient, scipy Voigt function and explicit Gaussian convolution. Saved objectives/Jacobians/arrays, stopping flags, AICc selection and conditional adequacy gates are replayed.','',f'Maximum finite-difference relative L2 derivative error: {result["max_finite_difference_relative_l2_error"]}.','', '|Target|Selected fit|χ² at21|χ²41−21|Maximum flux difference / error|','|---|---|---:|---:|---:|']
 for q in QUADRATURE:lines.append(f'|{q["target"]}|{q["name"]}|{q["chi2_21"]:.6f}|{q["delta_chi2_41_minus21"]:.6g}|{q["max_profile_change_in_error_units"]:.6g}|')
 lines+=['', 'Quadrature checks reuse fixed fitted parameters; they are not new optimizations. Successful ftol/xtol stops do not imply gtol stationarity or a global optimum. The separate stationarity/projection audit records first-order values and local conditional nuisance-projected errors. These are not calibrated discovery significances or alpha/time-law measurements.','', 'Sparse-profile provenance: the original generic driver saved `initial_parameters` before inserting the fixed2382 value. The actual full initial vector is reconstructed by setting that one entry to `fixed_2382_shift_km_s`; active-coordinate starts and fit results are unchanged. This audit preserves source/results and records the caveat.']
 lines+=['', 'The actual per-count selector uses the lowest chi-square successful stop when available, otherwise the lowest incomplete endpoint. The latter fallback is recorded explicitly in `validation/per_count_selection_audit.json`; every lower unfinished endpoint is retained. Later bounded recovery fits are audited as additional products, without rewriting original selection history.','', '|Target|Stage|Count|Selected endpoint|Lowest unfinished endpoint|Selected minus unfinished χ²|','|---|---|---:|---|---|---:|']
 for row in COUNT_SELECTION:
  if row['lower_incomplete_available']:lines.append(f'|{row["target"]}|{row["stage"]}|{row["ncomp"]}|{row["selected_name"]}|{row["lowest_incomplete_name"]}|{row["selected_minus_lowest_incomplete_chi2"]:.6f}|')
 for warning in WARNINGS:
  if warning['kind']=='lower_incomplete_conventional_endpoint':lines+=['',f'**Optimization warning for {warning["target"]}:** unfinished `{warning["name"]}` evaluated at the same OS21 quadrature has χ²={warning["endpoint_chi2_at21"]:.9f}, {warning["below_selected_null"]:.6f} below the selected H0 and {warning["below_selected_alternative"]:.6f} below selected H1. The reported pair is a conditional basin comparison, not evidence of improvement over the best available conventional endpoint.']
 if result['failures']:lines+=['','Failed checks:','']+[f'- {x["name"]}: {x["detail"]}' for x in result['failures']]
 (BASE/'reports/archive_complete_validation.md').write_text('\n'.join(lines)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['checks','quadrature','warnings','selected']},indent=2),flush=True)

if __name__=='__main__':main()
