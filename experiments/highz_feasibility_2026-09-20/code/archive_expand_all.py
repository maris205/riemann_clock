#!/usr/bin/env python3
"""Retrieve all unique sightlines passing the frozen primary metadata screen."""
import concurrent.futures as cf
import json
from datetime import datetime, timezone
from archive_fetch import ROOT, OUT, fetch_one, sha

def main():
    screen=json.loads((OUT/'catalogue_screen.json').read_text())
    candidates=[r for r in screen['all_absorbers'] if r['metadata_eligible']]
    unique={}
    for r in candidates:
        unique.setdefault(r['target'],{k:r[k] for k in ['target','z_abs','rank_score','n_strong_eligible']})
    protocol={'frozen_utc':datetime.now(timezone.utc).isoformat(),
       'rationale':'User requested all presently feasible archive analyses. Extend downloaded quality audit to every absorber passing original metadata criteria, before any line-offset measurement; preserve initial top5 and secondarytop3 selection history.',
       'n_absorbers':len(candidates),'n_unique_sightlines':len(unique),
       'source_screen_sha256':sha((OUT/'catalogue_screen.json').read_bytes()),
       'absorbers':[{'target':r['target'],'z_abs':r['z_abs']} for r in candidates]}
    (OUT/'complete_scope_policy.json').write_text(json.dumps(protocol,indent=2)+'\n')
    resources=json.loads((ROOT/'data/raw/catalog_portal_metadata.json').read_text())['result']['resources']
    lookup={r['name']:r for r in resources};results=[];failures=[]
    with cf.ThreadPoolExecutor(max_workers=5) as exe:
        futures={exe.submit(fetch_one,s,lookup[s['target']+'_Final_Spectrum.tar.gz']):s['target'] for s in unique.values()}
        for f in cf.as_completed(futures):
            name=futures[f]
            try:
                r=f.result();results.append(r)
                print(name,'ok',r['archive_bytes'],flush=True)
            except Exception as e:
                failures.append({'target':name,'error':repr(e)});print(name,'FAILED',repr(e),flush=True)
    result={'created_utc':datetime.now(timezone.utc).isoformat(),
        'scope_policy_sha256':sha((OUT/'complete_scope_policy.json').read_bytes()),
        'n_absorbers':len(candidates),'n_unique_sightlines':len(unique),'sources':sorted(results,key=lambda r:r['target']),
        'failures':failures,'total_archive_bytes':sum(r['archive_bytes'] for r in results),
        'license':'Spectrum portal CC BY-SA, version unspecified; metadata/code CC BY4.0. Cite Murphy et al.2019 DOI10.1093/mnras/sty2834.'}
    (OUT/'complete_source_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='sources'},indent=2))
if __name__=='__main__':main()
