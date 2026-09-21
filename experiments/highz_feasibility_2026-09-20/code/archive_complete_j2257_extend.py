#!/usr/bin/env python3
"""One bounded null-model count extension, justified before any free offsets."""
import concurrent.futures as cf,json
from datetime import datetime,timezone
from archive_complete_model import OUT,read_json,fit,run_alternative,sha
FOLDER=OUT/'J225719-100104'
def job(args):return fit(**args)
def main():
 config=read_json(FOLDER/'selection_frozen.json');oldsummary=read_json(FOLDER/'ordinary_model_summary.json')
 assert not oldsummary['passes_conditional_adequacy_gate']
 assert not (FOLDER/'selected_pair.json').exists(), 'Extension must precede relative-shift model selection'
 protocol={'frozen_utc':datetime.now(timezone.utc).isoformat(),
  'rationale':'Original null AICc falls continuously through22components (bestH0chi683.75/426), beforeanyshiftfit. Extend once to26/30 to check finite conventional-gas underfitting; stop at30 regardless outcome.',
  'counts':[26,30],'seeds':[0,1,2],'initial_max_nfev':350,'unchanged_configuration_sha256':sha(FOLDER/'selection_frozen.json'),
  'original_null_summary_sha256':sha(FOLDER/'ordinary_model_summary.json'),
  'selection':'Compare null AICc only; do not select by any relative displacement or alternative gain.'}
 (FOLDER/'count_extension_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
 jobs=[dict(config=config,ncomp=n,free=False,seed=s,name=f'n{n}_null_extension_s{s}_os9',oversample=9,max_nfev=350) for n in[26,30] for s in[0,1,2]]
 with cf.ProcessPoolExecutor(3) as pool:records=list(pool.map(job,jobs))
 original=[read_json(FOLDER/(r['name']+'.json')) for r in oldsummary['by_count']]
 selected=min([r for r in records+original if r['optimizer_success']],key=lambda r:r['aicc'])
 refined=fit(config,selected['ncomp'],False,selected['seed'],selected,name=f'n{selected["ncomp"]}_null_extension_refined_os21',oversample=21,max_nfev=600)
 ordinary={'target':config['target'],'selected_null':refined['name'],'ncomp':refined['ncomp'],
  'passes_conditional_adequacy_gate':bool(refined['optimizer_success'] and refined['chi2_per_ndf']<=1.5 and max(x['chi2_per_pixel'] for x in refined['per_line'])<=1.8),
  'extension_count_records':[{'name':r['name'],'ncomp':r['ncomp'],'chi2':r['chi2'],'aicc':r['aicc'],'optimizer_success':r['optimizer_success']} for r in records],
  'unchanged_original_summary':'ordinary_model_summary.json','extension_protocol':'count_extension_protocol.json'}
 (FOLDER/'ordinary_model_extension_summary.json').write_text(json.dumps(ordinary,indent=2)+'\n')
 if ordinary['passes_conditional_adequacy_gate']:run_alternative(config,refined,ordinary)
 else:print('Finite extension remains inadequate; no free offsets.',flush=True)
 print(json.dumps(ordinary,indent=2))
if __name__=='__main__':main()
