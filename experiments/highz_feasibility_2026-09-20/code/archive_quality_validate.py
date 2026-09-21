import csv,json,math,hashlib
from pathlib import Path
b=Path(__file__).resolve().parents[1]
j=json.loads((b/'results/archive_expansion/catalogue_screen.json').read_text()); p=j['policy']; checks=[]
def check(name,ok,detail=None):
 checks.append({'name':name,'passed':bool(ok),'detail':detail})
def close(a,c):return abs(a-c)<=1e-8*max(1,abs(c))
qs={x['Name_Adopt']:x for x in csv.DictReader((b/'data/raw/DR1_quasars_master.csv').open())}
ds=list(csv.DictReader((b/'data/raw/DR1_DLAs.csv').open()))
atomic={}
for line in (b/'data/raw/espresso_null/MM_VPFIT_2013-11-10_noiso.dat').read_text().splitlines():
 cols=line.split()
 if cols and cols[0]=='FeII':
  lam=float(cols[1]);atomic[str(int(lam))]=(lam,float(cols[2]))
rows=[]; total_lines=0
for d in ds:
 q=qs[d['Name_Adopt']];zem=float(q['zem_Adopt']); cover=[tuple(map(float,s.split('-'))) for s in q['WavCoverage'].split(',') if s];cnrs=list(map(float,q['CNR'].split(',')))
 for zstr in d['DLAzabs'].split(','):
  z=float(zstr); sr=next(x for x in j['all_absorbers'] if x['target']==d['Name_Adopt'] and x['z_abs']==z); eligible=[]
  for sl in sr['lines']:
   total_lines+=1; name=f'{d["Name_Adopt"]}:{z}:{sl["line"]}';lam,f=atomic[sl['line']]; wl=lam*(1+z);lo=wl*math.exp(-p['half_window_km_s']/299792.458);hi=wl*math.exp(p['half_window_km_s']/299792.458);nearest=min(range(5),key=lambda k:abs([3500,4500,5500,6500,7500][k]-wl));cnr=cnrs[nearest]; reasons=[]
   if not any(x<=lo and hi<=y for x,y in cover): reasons.append('catalogue_coverage_gap')
   if not(p['observed_wavelength_range_AA'][0]<=lo and hi<=p['observed_wavelength_range_AA'][1]):reasons.append('outside_fixed_wavelength_range')
   if any(max(lo,x)<min(hi,y) for x,y in p['conservative_telluric_exclusion_AA']):reasons.append('conservative_telluric_band')
   if lo<=1215.6701*(1+zem)*math.exp(p['red_of_quasar_Lyalpha_km_s']/299792.458):reasons.append('forest_or_Lyalpha_proximity')
   if cnr<p['metadata_CNR_per_2p5_km_s_min']:reasons.append('low_catalogue_CNR_proxy')
   check(name+' laboratory',close(sl['rest_AA'],lam) and close(sl['f'],f));check(name+' wavelength/window',close(sl['observed_AA'],wl) and close(sl['window_AA'][0],lo) and close(sl['window_AA'][1],hi));check(name+' CNR',sl['CNR_proxy']==cnr);check(name+' gates',sorted(sl['rejections'])==sorted(reasons) and sl['eligible']==(len(reasons)==0))
   if not reasons and f>=.01:eligible.append(cnr)
  reasons=[]
  if len(eligible)<p['minimum_strong_lines_f_ge_0p01']:reasons.append('fewer_than_3_eligible_strong_lines')
  if 299792.458*math.log((1+zem)/(1+z))<p['quasar_proximity_km_s']:reasons.append('absorber_quasar_proximity')
  if d['Name_Adopt'] in p['previous_targets_excluded']:reasons.append('already_analyzed_target')
  score=sorted(eligible,reverse=True)[2] if len(eligible)>=3 else 0
  check(f'{d["Name_Adopt"]}:{z}:absorber gate',reasons==sr['rejections'] and sr['metadata_eligible']==(not reasons));check(f'{d["Name_Adopt"]}:{z}:rank',len(eligible)==sr['n_strong_eligible'] and score==sr['rank_score'])
  if not reasons:rows.append({'target':d['Name_Adopt'],'z_abs':z,'rank_score':score,'n_strong_eligible':len(eligible)})
ranked=sorted(rows,key=lambda r:(-r['rank_score'],-r['n_strong_eligible'],r['target'],r['z_abs']));selected=[]; seen=set()
for r in ranked:
 if r['target'] not in seen: selected.append(r);seen.add(r['target'])
 if len(selected)==5:break
check('all catalogue rankings',ranked==j['ranked_absorbers']);check('five distinct sightlines',selected==j['selected']);check('155 absorbers',sum(len(x['DLAzabs'].split(',')) for x in ds)==155==j['DLA_absorbers']);check('36 metadata eligible',len(rows)==36==j['metadata_eligible_absorbers']);check('132 DLA sightlines',len(ds)==132==j['DLA_sightlines']);check('467 final spectra',sum(q['Spec_status'] in ['0','1','1,2'] for q in qs.values())==467==j['final_spectra_positive_exposure_count'])
for name,sha in j['source_hashes'].items():check('source hash '+name,hashlib.sha256((b/name).read_bytes()).hexdigest()==sha)
out={'scope':'Independent catalogue reconstruction; no screening-code imports','passed':all(c['passed'] for c in checks),'n_checks':len(checks),'n_absorbers':155,'n_line_windows':total_lines,'n_eligible':len(rows),'selected':selected,'failures':[c for c in checks if not c['passed']],'checks':checks}
(b/'results/archive_expansion/independent_catalogue_validation.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
