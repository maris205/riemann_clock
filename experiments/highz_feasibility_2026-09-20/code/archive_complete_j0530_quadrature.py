"""Fixed-parameter pixel-quadrature check; not an additional optimization."""
from pathlib import Path
import json
import numpy as np
from archive_complete_model import Model
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/archive_complete/J053007-250329'
s=json.loads((OUT/'selected_pair.json').read_text());out={'method':'Replay identical saved parameters at higher pixel-integration subdivisions; no reoptimization.','cases':{}}
for key in ['null','alternative']:
 r=json.loads((OUT/(s[key]+'.json')).read_text());values={}
 for os in [9,21,33,55]:
  m=Model(r['configuration'],r['ncomp'],r['free_shifts'],oversample=os);res,_,pro=m.evaluate(np.array(r['parameters']));values[str(os)]={'chi2':float(res@res),'profiles':pro}
 reference=values['55']
 out['cases'][key]={'name':r['name'],'oversamples':{o:{'chi2':q['chi2'],'chi2_minus_os55':q['chi2']-reference['chi2'],'max_abs_flux_vs_os55':float(max(np.max(np.abs(a-b)) for a,b in zip(q['profiles'],reference['profiles'])))} for o,q in values.items()}}
for os in [9,21,33,55]:
 q0=out['cases']['null']['oversamples'][str(os)]['chi2'];q1=out['cases']['alternative']['oversamples'][str(os)]['chi2'];out.setdefault('delta_chi2_by_os',{})[str(os)]=q0-q1
(OUT/'quadrature_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
