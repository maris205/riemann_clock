#!/usr/bin/env python3
"""Predeclared fixed-architecture omission of sky-proxy-flagged weak 1611."""
from pathlib import Path
import json
from datetime import datetime,timezone
from archive_complete_model import OUT,read_json,fit,local_pair_uncertainty
TARGET='J064326-504112'
FOLDER=OUT/TARGET

def run():
 config=read_json(FOLDER/'config.json');ordinary=read_json(FOLDER/'ordinary_model_summary.json')
 if not ordinary['passes_conditional_adequacy_gate']:
  out={'target':TARGET,'run':False,'reason':'Primary conventional model did not pass its predeclared conditional adequacy gate; no relative-shift sensitivity is interpreted.'}
  (FOLDER/'omit1611_summary.json').write_text(json.dumps(out,indent=2)+'\n');return out
 primary=read_json(FOLDER/'selected_pair.json');start=read_json(FOLDER/(primary['null']+'.json'));config=dict(config);config['lines']=[1608,2374,2382]
 config['sensitivity']='Predeclared omission of sky-proxy-flagged weak 1611 at fixed primary architecture, window, and retained native pixels. No component reselection or residual clipping.'
 (FOLDER/'omit1611_config.json').write_text(json.dumps(config,indent=2)+'\n')
 n=primary['ncomp'];seed=start['seed'];null=fit(config,n,False,seed,start,name=f'omit1611_n{n}_null_os21',oversample=21,max_nfev=500)
 alt=fit(config,n,True,seed,null,name=f'omit1611_n{n}_shift_os21',oversample=21,max_nfev=500)
 cross=fit(config,n,False,seed,alt,name=f'omit1611_n{n}_null_cross_os21',oversample=21,max_nfev=500)
 best=min([null,cross],key=lambda r:r['chi2'])
 if best['chi2']<null['chi2']-1e-3:
  alt2=fit(config,n,True,seed,best,name=f'omit1611_n{n}_shift_cross_os21',oversample=21,max_nfev=500);alt=min([alt,alt2],key=lambda r:r['chi2'])
 out={'target':TARGET,'run':True,'ncomp':n,'ndata':best['ndata'],'null':best['name'],'alternative':alt['name'],'chi2_null':best['chi2'],'chi2_alternative':alt['chi2'],'nominal_ndf_null':best['nominal_ndf'],'delta_chi2':best['chi2']-alt['chi2'],'extra_shifts':alt['npar']-best['npar'],'both_optimizer_success':bool(best['optimizer_success'] and alt['optimizer_success']),'local_pair_uncertainty':local_pair_uncertainty(alt),'null_per_line':best['per_line'],'null_stationarity':best['stationarity'],'alternative_stationarity':alt['stationarity'],'interpretation':'Fixed-primary-architecture sky-line sensitivity; losing weak1611 saturation information is confounded with omitting possible sky contamination. Diagonal expected errors; no precision calibration, alpha, or time-law measurement.'}
 (FOLDER/'omit1611_summary.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':print(json.dumps(run(),indent=2))
