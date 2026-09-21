#!/usr/bin/env python3
"""Bounded diagnostic: can shifts repair the failed four-line primary model?

Protocol amendment authorized after the quality-gated campaign. This diagnostic
cannot promote J0643 to a precision sample or substitute for calibrated tests.
"""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib
import archive_complete_model as engine
ROOT=Path(__file__).resolve().parents[1]
TARGET='J064326-504112';PRIMARY=ROOT/'results/archive_complete'/TARGET
OUT=PRIMARY/'failure_pair'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 config=engine.read_json(PRIMARY/'config.json');ordinary=engine.read_json(PRIMARY/'ordinary_model_summary.json');null=engine.read_json(PRIMARY/(ordinary['selected_null']+'.json'))
 assert config['lines']==[1608,1611,2374,2382] and null['ndata']==436 and null['ncomp']==16
 existing=[p.name for p in (OUT/TARGET).glob('*.json')] if (OUT/TARGET).exists() else []
 if existing:raise RuntimeError('Refusing to overwrite an existing failure-pair run')
 amendment={'declared_utc':datetime.now(timezone.utc).isoformat(),'authorization':'Root-authorized bounded failure-case diagnostic after the primary quality screen. Assess whether excluding all inadequate nulls would hide sensitivity to real frequency changes.','original_primary':null['name'],'original_primary_sha256':sha(PRIMARY/(null['name']+'.json')),'original_primary_arrays_sha256':sha(PRIMARY/(null['name']+'.npz')),'original_config_sha256':sha(PRIMARY/'config.json'),'configuration':config,'ndata':436,'ncomp':16,'expected_error_row':2,'oversample':21,'sequence':['One H1 from best saved primary H0; max600 evaluations.','One H0 cross-start from that H1; max500 evaluations.','Only if H0 improves by more than0.01 in chi2: one H1 cross-start from improved H0, max500 evaluations.'],'shift_anchor':2374,'free_shift_keys':[1608,1611,2382],'shift_bounds_km_s':[-1,1],'predeclared_interpretation':'Report shifts, model adequacy, active bounds, first-order diagnostics and the objective improvement even if failed. Initial gate was a quality screen, not a population null test. No discovery significance, precision displacement, alpha estimate or cosmic-time point. Primary remains gated out regardless of delta.','original_files_protected':True,'runner_sha256':sha(__file__),'engine_sha256':sha(engine.__file__)}
 write(OUT/'protocol_amendment.json',amendment)
 # The shared engine appends the unchanged physical target name. This output
 # root isolates every new checkpoint, record and array beneath failure_pair/.
 engine.OUT=OUT
 alt=engine.fit(config,16,True,null['seed'],null,name='failure_n16_shift_from_primary_os21',oversample=21,max_nfev=600)
 cross=engine.fit(config,16,False,null['seed'],alt,name='failure_n16_null_cross_os21',oversample=21,max_nfev=500)
 bestnull=min([null,cross],key=lambda r:r['chi2']);alts=[alt]
 if bestnull['chi2']<null['chi2']-.01:
  alts.append(engine.fit(config,16,True,bestnull['seed'],bestnull,name='failure_n16_shift_cross_os21',oversample=21,max_nfev=500))
 bestalt=min(alts,key=lambda r:r['chi2'])
 def gate(r):return bool(r['optimizer_success'] and r['chi2_per_ndf']<=1.5 and max(x['chi2_per_pixel'] for x in r['per_line'])<=1.8)
 summary={'target':TARGET,'classification':'post-quality-screen failed-primary explanatory comparison; not a precision observation','primary_remains_gated_out':True,'common_observable_status':'not_estimated_primary_due_model_failure','original_primary_null':null['name'],'original_primary_null_chi2':null['chi2'],'original_primary_null_relative_path':'../'+null['name']+'.json','null':bestnull['name'],'null_relative_path':str(Path(TARGET)/(bestnull['name']+'.json')) if bestnull is cross else '../'+null['name']+'.json','alternative':bestalt['name'],'alternative_relative_path':str(Path(TARGET)/(bestalt['name']+'.json')),'ndata':436,'ncomp':16,'npar_null':bestnull['npar'],'npar_alternative':bestalt['npar'],'chi2_null':bestnull['chi2'],'chi2_alternative':bestalt['chi2'],'nominal_ndf_null':bestnull['nominal_ndf'],'nominal_ndf_alternative':bestalt['nominal_ndf'],'delta_chi2':bestnull['chi2']-bestalt['chi2'],'delta_from_original_null':null['chi2']-bestalt['chi2'],'extra_shift_parameters':3,'shifts_m_s':bestalt['shifts_m_s'],'both_optimizer_success':bool(bestnull['optimizer_success'] and bestalt['optimizer_success']),'null_passes_original_conditional_gate':gate(bestnull),'alternative_passes_same_residual_gate':gate(bestalt),'null_per_line':bestnull['per_line'],'alternative_per_line':bestalt['per_line'],'null_active_bounds':bestnull['active_bounds'],'alternative_active_bounds':bestalt['active_bounds'],'null_stationarity':bestnull['stationarity'],'alternative_stationarity':bestalt['stationarity'],'original_primary_hash_unchanged':sha(PRIMARY/(null['name']+'.json'))==amendment['original_primary_sha256'],'original_config_hash_unchanged':sha(PRIMARY/'config.json')==amendment['original_config_sha256'],'interpretation':amendment['predeclared_interpretation']}
 write(OUT/'summary.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
