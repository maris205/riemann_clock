#!/usr/bin/env python3
"""Pre-fit deterministic source, neighbor-opacity, uniqueness and Jacobian checks."""
from pathlib import Path
import json, numpy as np
from archive_complete_j2331_model import NeighborModel, CoreModel, engine
ROOT=engine.ROOT;F=ROOT/'results/archive_complete/J233156-090802'
cfg=engine.read_json(F/'selection.json');m=NeighborModel(cfg,6,True,9);bare=CoreModel(cfg,6,True,9)
p,lo,hi=m.initial(0);p=np.maximum(lo+1e-4,np.minimum(hi-1e-4,p));r,J,prof=m.evaluate(p)
checks={}
checks['three_likelihood_regions']=len(m.lines)==3 and [L['key'] for L in m.lines]==[1608,2374,2382]
indices=np.concatenate([L['source_indices'][L['good']] for L in m.lines]);checks['global_source_indices_unique']=len(indices)==len(np.unique(indices))==1176
checks['one_1608_plus_1611_opacity_list']=[len(L['atom']) for L in m.lines]==[8,4,4]
checks['primary_references_unchanged']=all(a['ref']==b['ref'] for a,b in zip(m.lines,bare.lines))
checks['pixel_arrays_unchanged']=all(np.array_equal(a[k],b[k]) for a,b in zip(m.lines,bare.lines) for k in ['v','wave','grid','flux','error','good','source_indices'])
checks['adapter_hash_matches']=cfg['neighbor_opacity_adapter_sha256']==engine.sha(ROOT/cfg['neighbor_opacity_adapter'])
rb,Jb,profb=bare.evaluate(p);checks['other_region_models_unchanged']=all(np.array_equal(a,b) for a,b in zip(prof[1:],profb[1:]))
checks['neighbor_opacity_effect_nonzero']=float(np.max(np.abs(prof[0]-profb[0])))>.001
checks['neighbor_opacity_does_not_add_parameters']=m.labels==bare.labels
jac=[]
for j,label in enumerate(m.labels):
 h=1e-5*max(1,abs(p[j]));pp=p.copy();pm=p.copy();pp[j]+=h;pm[j]-=h
 fd=(m.fun(pp)-m.fun(pm))/(2*h);err=np.linalg.norm(fd-J[:,j])/max(np.linalg.norm(fd),np.linalg.norm(J[:,j]),1e-12)
 jac.append({'label':label,'relative_norm_error':float(err)});checks['jac_'+label]=bool(err<2e-4)
report={'pass':all(checks.values()),'checks':checks,'n_checks':len(checks),'ndata':m.ndata,'unique_source_pixels':len(np.unique(indices)),'maximum_jacobian_relative_norm_error':max(x['relative_norm_error'] for x in jac),'jacobian_details':jac,'configuration_sha256':engine.sha(F/'selection.json'),'neighbor_adapter_sha256':engine.sha(ROOT/'code/archive_complete_j2331_model.py'),'shared_engine_sha256':engine.sha(ROOT/'code/archive_complete_model.py'),'neighbor_model_max_flux_effect_at_initial':float(np.max(np.abs(prof[0]-profb[0]))),'interpretation':'Implementation/source checks only; no gas fit adequacy or new-physics conclusion.'}
(F/'preflight_validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['checks','jacobian_details']},indent=2));assert report['pass']
