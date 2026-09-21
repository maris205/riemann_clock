#!/usr/bin/env python3
"""Reconstruct the published manual-mask subset from the original UPL.
The extraction uses these rectangles conservatively on native-bin overlaps;
this parser does not claim to replay UVES_popler's redispersed flux actions.
"""
from pathlib import Path
import argparse,collections,hashlib,json,re
ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'data/raw/exposures'
REGIONS={'2374':(5105.59,5109.05),'2382':(5123.39,5127.35),'2600':(5590.81,5595.12)}

def reconstruct():
 text=(ROOT/'data/raw/espresso_null/hes0515m4414.upl').read_text();actions=[];types=collections.Counter()
 for block in re.split(r'(?m)^(?=\d+_ACTN)',text)[1:]:
  number,kind=map(int,re.search(r'(\d+)_ACTN = (\d+)',block).groups());types[kind]+=1
  if kind not in (1,2,3,4):continue
  coords=list(map(float,re.search(r'_CORD =\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)',block).groups()))
  labels=[k for k,(lo,hi) in REGIONS.items() if coords[0]<=hi and coords[1]>=lo]
  if not labels:continue
  row=re.search(r'_SPEC = (\d+)\s+\d+_ORDR = (\d+)',block)
  actions.append(dict(action_number=number,action_type=kind,wavelength_lower_A=coords[0],wavelength_upper_A=coords[1],flux_lower=coords[2],flux_upper=coords[3],exposure_one_based=int(row[1]) if row else None,row_one_based=int(row[2]) if row else None,overlapping_lines=labels))
 artifact=dict(interpretation='These are published manual clipping rectangles on REDISPERSED trace fluxes; native wavelength-overlap masking is a conservative transfer, not exact replay. One-based exposure and row correspond directly to UPL files and FITS rows + 1. No unclip actions in entire UPL.',source='../espresso_null/hes0515m4414.upl',actions=actions)
 return artifact,types

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--write-masks',action='store_true',help='Recreate the deterministic mask JSON, otherwise verify existing artifact without editing it.');args=parser.parse_args();artifact,types=reconstruct();path=RAW/'provenance_upl_line_mask_actions.json'
 if args.write_masks:path.write_text(json.dumps(artifact,indent=2))
 same=json.loads(path.read_text())==artifact
 out=dict(pass_=same,source_upl_sha256=hashlib.sha256((ROOT/'data/raw/espresso_null/hes0515m4414.upl').read_bytes()).hexdigest(),mask_artifact_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),full_upl_action_types=dict(types),selected_action_count=len(artifact['actions']),selected_action_counts_by_line=dict(collections.Counter(k for a in artifact['actions'] for k in a['overlapping_lines'])),has_any_unclip=bool(types[3] or types[4]),wavelength_frame='Imported VAC_BARY preserved by the official ESPRESSO reader and FTESPR branch of UVES_wpol; no additional heliocentric or air/vacuum conversion.',source_code_manifest='provenance_sources.json')
 (RAW/'provenance_mask_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));raise SystemExit(0 if same else 1)
