#!/usr/bin/env python3
"""Independent raw-FITS extraction and secondary-screen verification."""
import io, json, math, hashlib, sys, tarfile
from pathlib import Path
import numpy as np
B=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(B/'.deps'))
from astropy.io import fits
O=B/'results/archive_expansion';checks=[]
def check(n,v,d=None):checks.append({'name':n,'passed':bool(v),'detail':d})
def sha(x):return hashlib.sha256(x).hexdigest()
screen=json.loads((O/'catalogue_screen.json').read_text());strict=json.loads((O/'strict_catalogue_screen.json').read_text())
expected=[]
for r in screen['all_absorbers']:
 good=[]; warned=[]
 for l in r['lines']:
  bounds=l['window_AA'];pad=strict['heliocentric_padding_km_s'];lo=bounds[0]*math.exp(-pad/299792.458);hi=bounds[1]*math.exp(pad/299792.458)
  overlap=any(max(lo,b0)<min(hi,b1) for b0,b1 in strict['water_bands_AA'])
  if overlap:warned.append(l['line'])
  if l['eligible'] and l['f']>=.01 and not overlap:good.append(l)
 expected.append({'target':r['target'],'z_abs':r['z_abs'],'strict_eligible':r['metadata_eligible'] and len(good)>=3,'n_strong_eligible':len(good),'strong_lines':[x['line'] for x in good],'broad_water_warning_lines':warned,'rank_score':sorted([x['CNR_proxy'] for x in good],reverse=True)[2] if len(good)>=3 else 0})
for actual, exp in zip(strict['all_absorbers'],expected):check(f"strict {exp['target']}:{exp['z_abs']}",actual==exp)
rank=sorted([r for r in expected if r['strict_eligible']],key=lambda r:(-r['rank_score'],-r['n_strong_eligible'],r['target'],r['z_abs']))
check('strict rank',rank==strict['ranked_eligible']);check('strict eligible count',len(rank)==7==strict['eligible_count'])
seen={r['target'] for r in screen['selected']};chosen=[]
for r in rank:
 if r['target'] not in seen:chosen.append(r);seen.add(r['target'])
 if len(chosen)==3:break
check('strict additional selection',chosen==strict['selected_additional'])
for key,path in [('catalogue_screen_sha256',O/'catalogue_screen.json'),('source_paper_sha256',B/'data/raw/paper_accepted_2018-10-16.tex'),('atomic_data_sha256',B/'data/raw/espresso_null/MM_VPFIT_2013-11-10_noiso.dat')]:check('strict provenance '+key,sha(path.read_bytes())==strict[key])
main=json.loads((O/'source_manifest.json').read_text()); sec=json.loads((O/'strict_source_manifest.json').read_text()); sources=main['sources']+sec['sources']
check('initial frozen-policy hash',sha((O/'selection_policy.json').read_bytes())==main['selection_policy_sha256']);check('initial screen hash',sha((O/'catalogue_screen.json').read_bytes())==main['catalogue_screen_sha256']);check('strict manifest hash',sha((O/'strict_catalogue_screen.json').read_bytes())==sec['strict_screen_sha256']);check('initial eight distinct targets',len({r['target'] for r in sources})==8==len(sources))
if (O/'complete_source_manifest.json').exists():
 complete=json.loads((O/'complete_source_manifest.json').read_text());scope=json.loads((O/'complete_scope_policy.json').read_text());sources=complete['sources']
 check('complete scope hash',sha((O/'complete_scope_policy.json').read_bytes())==complete['scope_policy_sha256']);check('complete source screen hash',sha((O/'catalogue_screen.json').read_bytes())==scope['source_screen_sha256'])
 wanted=sorted((r['target'],r['z_abs']) for r in screen['all_absorbers'] if r['metadata_eligible']);actual=sorted((r['target'],r['z_abs']) for r in scope['absorbers'])
 check('complete all36absorbers',wanted==actual and len(actual)==36==complete['n_absorbers']);check('complete all32sightlines',{r['target'] for r in sources}=={r[0] for r in wanted} and len(sources)==32==complete['n_unique_sightlines']);check('complete zero download failures',len(complete['failures'])==0);check('complete total archive bytes',sum(r['archive_bytes'] for r in sources)==complete['total_archive_bytes'])
for r in sources:
 n=r['target'];path=B/r['archive_file'];blob=path.read_bytes();check(n+' archive hash',sha(blob)==r['archive_sha256']);check(n+' archive length',len(blob)==r['archive_bytes'])
 with tarfile.open(fileobj=io.BytesIO(blob),mode='r:gz') as t:
  member=t.getmember(f'{n}/{n}.fits');fb=t.extractfile(member).read()
 check(n+' FITS hash',sha(fb)==r['fits_sha256'])
 with fits.open(io.BytesIO(fb)) as hdus:
  a=hdus[0].data;h=hdus[0].header;check(n+' 9-row schema',a.shape[0]==9);check(n+' log scale',h['DC-FLAG']==1 and h['CTYPE1']=='LOGLIN');check(n+' FITS dimension',h['NAXIS1']==a.shape[1]==r['pixels'])
  wl=np.power(10.,h['CRVAL1'])*np.power(10.,(np.arange(a.shape[1])+1-h['CRPIX1'])*h['CD1_1']);npz=np.load(B/'data/processed/archive_expansion'/f'{n}.npz')
  check(n+' wavelengths',np.allclose(wl,npz['wavelength_AA'],rtol=3e-15,atol=0));check(n+' monotonic',np.all(np.diff(wl)>0));check(n+' endpoints',np.allclose([wl[0],wl[-1]],r['wavelength_AA_range'],rtol=3e-15,atol=0))
  for k,row in [('flux',0),('error',1),('expected_fluctuation',2),('continuum',3),('status',4),('n_contrib_before',5),('n_contrib_after',6),('chi2_before',7),('chi2_after',8)]:check(n+' array '+k,np.array_equal(npz[k],a[row],equal_nan=True))
  valid=(a[4]==1)&np.isfinite(a[0])&np.isfinite(a[1])&(a[1]>0);check(n+' validity',np.array_equal(valid,npz['valid']));check(n+' valid count',int(valid.sum())==r['valid_pixels']);check(n+' exposure count',h['UP_NEXP']==r['n_exposures']);check(n+' dispersion',h['UP_DISP']==r['native_dispersion_km_s']);check(n+' elapsed exposure',h['UP_TEXP']==r['exposure_s'])
  dates=[str(x).strip() for x in hdus[3].data['UTDate']];check(n+' observation dates',[min(dates),max(dates)]==r['observation_date_range'])
result={'scope':'Independent secondary-screen arithmetic and arrays reconstructed directly from original FITS tar members; no analysis-code imports','passed':all(x['passed'] for x in checks),'n_checks':len(checks),'n_spectra':len(sources),'strict_eligible':len(rank),'failures':[x for x in checks if not x['passed']],'checks':checks}
(O/'independent_data_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
