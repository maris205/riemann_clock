#!/usr/bin/env python3
"""Bounded recovery of a materially lower capped ordinary-model endpoint.

This preserves the higher successfully-stopped pair as historical outputs. It
does not discard the lower unfinished null in favor of an apparent shift gain.
"""
import json, shutil
from datetime import datetime,timezone
import numpy as np
from archive_complete_model import OUT,read_json,fit,local_pair_uncertainty,sha
FOLDER=OUT/'J225719-100104'
def main():
 config=read_json(FOLDER/'selection_frozen.json')
 low=read_json(FOLDER/'n30_null_extension_s2_os9.json')
 before=read_json(FOLDER/'selected_pair.json')
 protocol={'frozen_utc':datetime.now(timezone.utc).isoformat(),
  'reason':'A capped same-architecture null endpoint has chi2=422.187928, below both earlier successful null463.718227 and alternative454.225607. Successful stopping must not discard this lower feasible endpoint.',
  'source_low_endpoint':low['name'],'source_sha256':sha(FOLDER/(low['name']+'.json')),
  'budget':'One600-evaluation null continuation atOS21, oneH1fromrecoverednull, oneH0cross, atmostoneH1cross, sparse5pointprofile300each. Noadditionalcountsorwindows.',
  'previous_pair':before}
 (FOLDER/'lower_endpoint_recovery_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
 for name in ['selected_pair.json','pair_profile.json']:
  path=FOLDER/name
  if path.exists():shutil.copyfile(path,FOLDER/(path.stem+'_before_lower_endpoint_recovery.json'))
 null=fit(config,30,False,2,low,name='n30_null_recovered_os21',oversample=21,max_nfev=600)
 adequate=bool(null['optimizer_success'] and null['chi2_per_ndf']<=1.5 and max(r['chi2_per_pixel'] for r in null['per_line'])<=1.8)
 if not adequate:
  result={'status':'lower_endpoint_continuation_not_adequate_or_incomplete','record':null['name'],
    'previous_pair_superseded_for_inference':True,'relative_shift_estimate_not_certified':True}
  (FOLDER/'recovery_outcome.json').write_text(json.dumps(result,indent=2)+'\n');return
 alt=fit(config,30,True,2,null,name='n30_shift_recovered_os21',oversample=21,max_nfev=600)
 cross=fit(config,30,False,2,alt,name='n30_null_recovered_cross_os21',oversample=21,max_nfev=500)
 null=min([null,cross],key=lambda r:r['chi2'])
 if null['chi2']<read_json(FOLDER/'n30_null_recovered_os21.json')['chi2']-1e-3:
  ac=fit(config,30,True,2,null,name='n30_shift_recovered_cross_os21',oversample=21,max_nfev=500)
  alt=min([alt,ac],key=lambda r:r['chi2'])
 local=local_pair_uncertainty(alt)
 summary={'target':config['target'],'z_abs':config['z'],'lines':config['lines'],'window':config['window'],
  'ncomp':30,'ndata':null['ndata'],'null':null['name'],'alternative':alt['name'],
  'chi2_null':null['chi2'],'chi2_alternative':alt['chi2'],'nominal_ndf_null':null['nominal_ndf'],
  'delta_chi2':null['chi2']-alt['chi2'],'extra_shifts':alt['npar']-null['npar'],
  'both_optimizer_success':bool(null['optimizer_success'] and alt['optimizer_success']),
  'ordinary_gate':adequate,'local_pair_uncertainty':local,
  'null_stationarity':null['stationarity'],'alternative_stationarity':alt['stationarity'],
  'recovery_protocol':'lower_endpoint_recovery_protocol.json',
  'supersedes':'selected_pair_before_lower_endpoint_recovery.json',
  'interpretation':'Conditional exploration after retaining and continuing a materially lower cappednull. Finite multistarts do not certify globaloptimum; notalpha/time discovery.'}
 (FOLDER/'selected_pair.json').write_text(json.dumps(summary,indent=2)+'\n')
 sig=local['rank_sensitivity'][1]['conditional_sigma_m_s'];center=local['shift_m_s']/1000
 points=[]
 for i,x in enumerate(sorted(set(float(np.clip(center+k*sig/1000,-.999,.999)) for k in [-2,-1,0,1,2]))):
  p=fit(config,30,True,2,alt,name=f'n30_profile2382_recovered_{i}_os21',oversample=21,max_nfev=300,fixed_shift=x)
  points.append({'value_m_s':1000*x,'chi2':p['chi2'],'delta_from_free':p['chi2']-alt['chi2'],'optimizer_success':p['optimizer_success'],'name':p['name']})
 (FOLDER/'pair_profile.json').write_text(json.dumps({'free_name':alt['name'],'free_chi2':alt['chi2'],'grid':points,
  'supersedes':'pair_profile_before_lower_endpoint_recovery.json',
  'interpretation':'Sparse nuisance-profile diagnostic. Lowerprofilepoint would be an optimization warning, not negativeconfidence statistic.'},indent=2)+'\n')
 (FOLDER/'recovery_outcome.json').write_text(json.dumps({'status':'recovered_pair_saved','selected_pair':'selected_pair.json','ordinary_success':null['optimizer_success'],'alternative_success':alt['optimizer_success']},indent=2)+'\n')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
