#!/usr/bin/env python3
"""Fetch all 17 public S2D exposure products from a pinned source commit."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,hashlib,time,requests
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/raw/exposures';OUT.mkdir(parents=True,exist_ok=True)
TREE=json.loads((ROOT/'data/raw/espresso_github_tree.json').read_text()); COMMIT=TREE['sha']
BASE=f'https://raw.githubusercontent.com/MTMurphy77/ESPRESSO_HE0515-4414/{COMMIT}/'
def fetch(item):
 path=item['path']; f=OUT/Path(path).name
 for attempt in range(4):
  try:
   if not f.exists() or f.stat().st_size!=item['size']:
    response=requests.get(BASE+path,timeout=(30,180));response.raise_for_status()
    if len(response.content)!=item['size']: raise ValueError('size mismatch')
    temp=f.with_suffix('.download');temp.write_bytes(response.content);temp.replace(f)
   data=f.read_bytes(); gitsha=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
   if gitsha!=item['sha']: raise ValueError('Git blob hash mismatch')
   row=dict(source_url=BASE+path,repository_path=path,commit=COMMIT,local_file=str(f.relative_to(ROOT)),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),git_blob_sha1=gitsha,git_blob_verified=True)
   print(f.name,'verified',flush=True);return row
  except Exception:
   if attempt==3: raise
   time.sleep(1+attempt)
items=[a for a in TREE['tree'] if a['path'].startswith('DRS_extracted_spectra/') and a['path'].endswith('.fits')]
with ThreadPoolExecutor(max_workers=4) as pool: rows=list(pool.map(fetch,items))
manifest=dict(repository='https://github.com/MTMurphy77/ESPRESSO_HE0515-4414',commit=COMMIT,source_data_doi='10.5281/zenodo.5512490',downloaded_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),exposure_count=len(rows),total_bytes=sum(x['bytes'] for x in rows),files=rows)
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
