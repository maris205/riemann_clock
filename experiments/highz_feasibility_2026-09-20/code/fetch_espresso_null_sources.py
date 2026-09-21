#!/usr/bin/env python3
"""Verify cached sources; fetch missing public files at their pinned URLs."""
from pathlib import Path
import json,hashlib,requests
ROOT=Path(__file__).resolve().parents[1];RAW=ROOT/'data/raw/espresso_null'
entries=json.loads((RAW/'manifest.json').read_text())['files']
spectrum=json.loads((RAW/'spectrum_manifest.json').read_text());spectrum['local_path']='hes0515m4414.fits';entries.append(spectrum)
for row in entries:
 p=RAW/Path(row['local_path']).name
 if not p.exists():
  r=requests.get(row['url'],timeout=180);r.raise_for_status();p.write_bytes(r.content)
 assert p.stat().st_size==row['bytes'],p
 assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],p
 print(p.name,'verified',p.stat().st_size)
