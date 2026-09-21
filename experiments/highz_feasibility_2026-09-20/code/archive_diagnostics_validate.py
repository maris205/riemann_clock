#!/usr/bin/env python3
"""Independent scalar/window reconstruction directly from archived FITS."""
import io,json,math,sys,tarfile,hashlib
from pathlib import Path
import numpy as np
B=Path(__file__).resolve().parents[1];O=B/'results/archive_expansion';D=B/'data/processed/archive_expansion';sys.path.insert(0,str(B/'.deps'))
from astropy.io import fits
C=299792.458;checks=[]
def ck(n,v,d=None):checks.append({'name':n,'passed':bool(v),'detail':d})
def same(a,b):
 if a is None:return b is None or not np.isfinite(b)
 if b is None:return False
 return bool(np.isclose(a,b,rtol=2e-10,atol=1e-10))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
audit=json.loads((O/'absorber_quality.json').read_text());summ=json.loads((O/'quality_summary.json').read_text());pol=json.loads((O/'actual_quality_policy.json').read_text());screen=json.loads((O/'catalogue_screen.json').read_text());sl={(a['target'],a['z_abs']):a for a in screen['all_absorbers']};strict=json.loads((O/'strict_catalogue_screen.json').read_text());cache={};linecount=0;gapchanges=[];readies=[]
for a in audit['results']:
 target=a['target'];z=a['z_abs'];base=sl[(target,z)];lab={l['line']:l for l in base['lines']}
 if target not in cache:
  with tarfile.open(B/'data/raw/archive_expansion'/f'{target}_Final_Spectrum.tar.gz','r:gz') as tar:
   bb=tar.extractfile(f'{target}/{target}.fits').read()
  with fits.open(io.BytesIO(bb)) as h:
   data=h[0].data.copy();head=h[0].header;table=h[1].data.copy()
  wave=10**(head['CRVAL1']+(np.arange(data.shape[1])+1-head['CRPIX1'])*head['CD1_1']);cache[target]=(data,head,wave,table)
 data,head,wave,table=cache[target];valid=(data[4]==1)&np.isfinite(data[0])&np.isfinite(data[1])&(data[1]>0);dv=C*math.log(10)*head['CD1_1']; primary=[];strong=[];weak=[];det=[];detgap=[]
 for r in a['line_rows']:
  linecount+=1;n=f'{target}:{z}:{r["line"]}';l=lab[r['line']];v=C*np.log(wave/(l['rest_AA']*(1+z)));ix=np.where(np.abs(v)<=250)[0];vv=v[ix];f=data[0,ix];e=data[1,ix];g=valid[ix];core=g&(np.abs(vv)<=100);side=g&(np.abs(vv)>=150);dl=wave[ix]*2*math.sinh(dv/(2*C))/(1+z)
  # Pixel-edge integration agrees with producer but derives dv from CD1_1.
  ew=np.sum((1-f[core])*dl[core]);er=np.linalg.norm(e[core]*dl[core]);ewf=np.sum((1-f[g])*dl[g]);erf=np.linalg.norm(e[g]*dl[g]);cnr=np.median(1/e[g]) if g.any() else None;frac=float(g.mean()) if len(g) else 0
  floor=np.maximum(.02,3*np.where(e>0,e,np.inf));tau=np.full(len(g),np.nan);tau[g]=-np.log(np.maximum(f[g],floor[g]));na=tau*3.768e14/(l['f']*l['rest_AA']);sat=core&(f<=np.maximum(.1,3*e));useful=core&(f>np.maximum(.1,3*e))&(f<.95)&(1-f>3*e)
  edges=[];run=0
  for good in g:
   if good:edges.append(run);run=0
   else:run+=1
  edges.append(run)
  lo=l['window_AA'][0]*math.exp(-30/C);hi=l['window_AA'][1]*math.exp(30/C);ov=lambda bands:any(max(lo,x)<min(hi,y) for x,y in bands);water=ov(strict['water_bands_AA']);sky=ov([(x-1,x+1) for x in [5578.,5891.,5897.]]);artifact=ov([(4713,4719),(4741,4747),(5237,5243),(5580,5800)])
  edge=g&(np.abs(vv)>200);edgesnr=np.sum((1-f[edge])*dl[edge])/np.linalg.norm(e[edge]*dl[edge]) if edge.any() else None
  nearest=min(table,key=lambda q:abs(q['Wavelength']-l['observed_AA']))
  expected={'actual_log_grid_spacing_km_s':dv,'native_pixels':len(ix),'valid_pixels':int(g.sum()),'valid_fraction':frac,'longest_invalid_run_pixels':max(edges,default=0),'native_CNR':cnr,'sideband_median_flux':np.median(f[side]) if side.any() else None,'core_EW_rest_AA':ew,'core_EW_diagonal_error_AA':er,'core_EW_diagonal_SNR':ew/er if er else None,'full_EW_rest_AA':ewf,'full_EW_diagonal_error_AA':erf,'core_to_full_EW_ratio':ew/ewf if ewf>0 else None,'core_low_flux_pixels':int(sat.sum()),'core_valid_pixels':int(core.sum()),'core_low_flux_fraction':sat.sum()/core.sum() if core.any() else None,'core_resolved_unsaturated_detection_pixels':int(useful.sum()),'core_AOD_integral_cm2':np.nansum(na[core])*dv,'full_window_edge_EW_diagonal_SNR':edgesnr,'nearest_nominal_R':nearest['NomResolPower'],'nearest_arc_R':nearest['ArcResolPower'],'median_contributing_extracted_pixels':np.median(data[6,ix][g]) if g.any() else None}
  for key,value in expected.items():ck(n+' '+key,same(r[key],value))
  ck(n+' warnings',water==r['broader_water_band_warning'] and sky==r['strong_sky_feature_warning'] and artifact==r['reported_instrument_artifact_warning'])
  reason=list(l['rejections']) if not l['eligible'] else []
  if frac<.98:reason.append('native_valid_fraction_below_0p98')
  if cnr is None or cnr<15:reason.append('native_CNR_below_15')
  if er<=0 or ew/er<5:reason.append('core_EW_diagonal_detection_proxy_below_5')
  keep=not reason;ks=keep and not water;ck(n+' rejectionlist',reason==r['rejections']);ck(n+' primary/strict gate',r['detected_primary_candidate']==keep and r['detected_strict_candidate']==ks)
  win=np.load(D/f'{target}_z{z:.3f}_FeII{r["line"]}_window.npz')
  for key,value in [('source_indices',ix),('velocity_km_s',vv),('wavelength_AA',wave[ix]),('flux',f),('error',e),('valid',g),('AOD_tau',tau),('AOD_Na_per_km_s',na),('core_mask',core),('uncertainty_floor_affected',g&(f<=floor))]:ck(n+' saved '+key,np.allclose(win[key],value,rtol=2e-12,atol=1e-12,equal_nan=True))
  if l['f']>=.01 and keep:primary.append(r['line'])
  if l['f']>=.01 and ks:
   strong.append(r['line']);grid=np.arange(-100,100.1,5);vp=vv[g];idx=ix[g];sn=(1-f[g])/e[g];score=np.interp(grid,vp,sn,left=0,right=0)>=3
   pos=np.searchsorted(vp,grid);allowed=(pos>0)&(pos<len(vp));inside=np.where(allowed)[0];allowed[inside]&=(idx[pos[inside]]-idx[pos[inside]-1])==1
   det.append(score);detgap.append(score&allowed)
   if r['line'] in ['2374','2586'] and useful.sum()>=5:weak.append(r['line'])
 common=int((np.sum(det,axis=0)>=3).sum()) if det else 0;cg=int((np.sum(detgap,axis=0)>=3).sum()) if detgap else 0
 key='common_5km_s_grid_points_with_3line_absorption' if 'common_5km_s_grid_points_with_3line_absorption' in a else 'common_5km_s_bins_with_3line_absorption'
 got=a[key];expect=cg if ('Both bracketing native pixels must be valid' in pol.get('common_grid_note','')) else common
 ck(f'{target}:{z} primary strong',primary==a['primary_detected_strong_lines']);ck(f'{target}:{z} strict strong',strong==a['strict_detected_strong_lines']);ck(f'{target}:{z} weakerline gate',weak==a['resolved_weak_line_candidates']);ck(f'{target}:{z} common sampling',got==expect,{'old_interpolation_count':common,'gap_aware_count':cg,'reported':got});ready=len(strong)>=3 and bool(weak) and expect>=3;ck(f'{target}:{z} readiness',ready==a['passes_computational_readiness'])
 if ready:readies.append((target,z))
 if common!=cg:gapchanges.append({'target':target,'z_abs':z,'old':common,'gap_aware':cg,'readiness_changed':(len(strong)>=3 and bool(weak) and common>=3)!=(len(strong)>=3 and bool(weak) and cg>=3)})
for path,h in summ['hashes'].items():ck('summary hash '+path,sha(B/path)==h)
ck('audit policy hash',sha(O/'actual_quality_policy.json')==audit['policy_sha256']);ck('all36audited',len(audit['results'])==36==summ['candidate_absorbers']);ck('all324windows',linecount==324==summ['actual_windows']);ck('readinesscount',len(readies)==summ['computational_readiness_count'])
result={'scope':'Independent native-FITS scalar, pixel-mask, AOD-floor and selection reconstruction; no producer-code imports','passed':all(x['passed'] for x in checks),'n_checks':len(checks),'n_spectra':len(cache),'n_absorbers':len(audit['results']),'n_windows':linecount,'readiness_candidates':readies,'interpolation_gap_sensitivity':gapchanges,'failures':[x for x in checks if not x['passed']],'checks':checks}
(O/'independent_diagnostics_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({**{k:v for k,v in result.items() if k not in ['checks','failures']},'n_failures':len(result['failures']),'first_failures':result['failures'][:8]},indent=2))
