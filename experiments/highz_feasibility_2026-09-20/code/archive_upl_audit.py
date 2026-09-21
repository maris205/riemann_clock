#!/usr/bin/env python3
"""Record explicit UPL VSHT fields; absence is not a complete calibration audit."""
import hashlib,json,re
from pathlib import Path
B=Path(__file__).resolve().parents[1];out=[]
for p in sorted((B/'data/raw/archive_expansion').glob('*.upl')):
 entries=[]
 for m in re.finditer(r'^\s*(\d+)_VSHT\s*=\s*([^\n]+)',p.read_text(),re.M):
  a=list(map(float,m[2].split()));entries.append({'input_index':int(m[1]),'values_as_recorded':a})
 out.append({'target':p.stem,'file':str(p.relative_to(B)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'n_input_spectral_records_with_VSHT':len(entries),'n_nonzero_first_or_second_values':sum(any(v!=0 for v in e['values_as_recorded'][:2]) for e in entries),'n_nonzero_slope_second_values':sum(e['values_as_recorded'][1]!=0 for e in entries),'entries':entries})
r={'purpose':'Record explicit archive processing settings, not calibrate or certify wavelength scales','note':'VSHT records refer to input spectral records/arms, not necessarily distinct science exposures. First and second values encode shift and slope; exact units and implementation documented in prior UVES_popler source audit. Zero values do not prove all calibration choices absent; nonzero values do not prove residual distortions eliminated.','n_targets':len(out),'n_with_nonzero_shift_or_slope':sum(t['n_nonzero_first_or_second_values']>0 for t in out),'targets':out}
(B/'results/archive_expansion/upl_wavelength_shift_audit.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='targets'},indent=2))
