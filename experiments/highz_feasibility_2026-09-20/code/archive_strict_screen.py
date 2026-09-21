#!/usr/bin/env python3
"""Secondary metadata screen, fixed before measuring any archive line offsets.

The primary pre-screen is preserved. This excludes all broader H2O bands listed
by the source paper, not just the strongest sub-bands. It is conservative and
may exclude genuinely useful lines; no telluric transmission model is supplied.
"""
import concurrent.futures as cf
import hashlib, json
from datetime import datetime, timezone
from archive_screen import ROOT, OUT, C
from archive_fetch import fetch_one, sha
import numpy as np

BANDS=[[6470.,6600.],[6830.,7450.],[7820.,8620.],[8780.,10000.]]
def main():
    screen=json.loads((OUT/'catalogue_screen.json').read_text())
    full=[]
    for r in screen['all_absorbers']:
        eligible=[]; warnings=[]
        for l in r['lines']:
            lo,hi=np.array(l['window_AA'])*np.exp(np.array([-30.,30.])/C)
            warned=any(lo<b and hi>a for a,b in BANDS)
            if warned: warnings.append(l['line'])
            if l['eligible'] and l['f']>=.01 and not warned: eligible.append(l)
        full.append({'target':r['target'],'z_abs':r['z_abs'],
            'strict_eligible':r['metadata_eligible'] and len(eligible)>=3,
            'n_strong_eligible':len(eligible),'strong_lines':[l['line'] for l in eligible],
            'broad_water_warning_lines':warnings,
            'rank_score':sorted([l['CNR_proxy'] for l in eligible],reverse=True)[2] if len(eligible)>=3 else 0})
    ranked=sorted([r for r in full if r['strict_eligible']],key=lambda r:(-r['rank_score'],-r['n_strong_eligible'],r['target'],r['z_abs']))
    old={r['target'] for r in screen['selected']}; selected=[]
    for r in ranked:
        if r['target'] not in old:
            selected.append(r); old.add(r['target'])
        if len(selected)==3: break
    result={'frozen_utc':datetime.now(timezone.utc).isoformat(),
        'rationale':'Source-paper broad H2O bands added after primary metadata screen, before inspecting any new line profiles or measuring offsets. Preserve initial five candidates, add top three distinct targets that avoid broad H2O bands.',
        'source':'Murphy et al.2019 DOI10.1093/mnras/sty2834, subsection Telluric features; locally pinned paper tex',
        'water_bands_AA':BANDS,'heliocentric_padding_km_s':30,
        'catalogue_screen_sha256':sha((OUT/'catalogue_screen.json').read_bytes()),
        'source_paper_sha256':sha((ROOT/'data/raw/paper_accepted_2018-10-16.tex').read_bytes()),
        'atomic_data_sha256':sha((ROOT/'data/raw/espresso_null/MM_VPFIT_2013-11-10_noiso.dat').read_bytes()),
        'atomic_basis':'Murphy/Berengut2014 terrestrial composite rest wavelengths and f-values; no individual isotope fitting in readiness diagnostics',
        'eligible_count':len(ranked),'selected_additional':selected,'ranked_eligible':ranked,'all_absorbers':full}
    (OUT/'strict_catalogue_screen.json').write_text(json.dumps(result,indent=2)+'\n')
    resources=json.loads((ROOT/'data/raw/catalog_portal_metadata.json').read_text())['result']['resources']
    lookup={r['name']:r for r in resources}
    sources=[]
    with cf.ThreadPoolExecutor(max_workers=3) as exe:
        futures=[exe.submit(fetch_one,s,lookup[s['target']+'_Final_Spectrum.tar.gz']) for s in selected]
        for f in cf.as_completed(futures):
            row=f.result(); sources.append(row); print(json.dumps({k:row[k] for k in ['target','z_abs','archive_bytes','n_exposures']}),flush=True)
    (OUT/'strict_source_manifest.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),'strict_screen_sha256':sha((OUT/'strict_catalogue_screen.json').read_bytes()),'sources':sources},indent=2)+'\n')
    print('strict eligible:',len(ranked),'additional:',len(sources))
if __name__=='__main__': main()
