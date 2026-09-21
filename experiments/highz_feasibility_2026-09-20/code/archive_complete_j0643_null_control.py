#!/usr/bin/env python3
"""Post-primary null-only diagnostic of a pre-fit sky-proxy-flagged weak line."""
from pathlib import Path
import json
from archive_complete_model import OUT,read_json,fit
FOLDER=OUT/'J064326-504112'
if __name__=='__main__':
 ordinary=read_json(FOLDER/'ordinary_model_summary.json');start=read_json(FOLDER/(ordinary['selected_null']+'.json'));config=read_json(FOLDER/'config.json');config=dict(config);config['lines']=[1608,2374,2382]
 config['sensitivity']='After failed four-line primary, fixed16 null-only diagnostic omitting weak1611. Sky-proxy warning was documented before primary. Same remaining native pixels/window/expected errors, start from primary gas; no architecture reselection or free shifts. This does not replace primary.'
 (FOLDER/'omit1611_null_only_config.json').write_text(json.dumps(config,indent=2)+'\n')
 r=fit(config,start['ncomp'],False,start['seed'],start,name=f'omit1611_n{start["ncomp"]}_null_only_os21',oversample=21,max_nfev=500)
 out={'name':r['name'],'target':r['target'],'ndata':r['ndata'],'ncomp':r['ncomp'],'chi2':r['chi2'],'nominal_ndf':r['nominal_ndf'],'optimizer_success':r['optimizer_success'],'per_line':r['per_line'],'passes_same_conditional_gate':bool(r['optimizer_success'] and r['chi2_per_ndf']<=1.5 and max(x['chi2_per_pixel'] for x in r['per_line'])<=1.8),'no_relative_shifts_fitted':True,'primary_still_inadequate':True,'interpretation':config['sensitivity']}
 (FOLDER/'omit1611_null_only_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
