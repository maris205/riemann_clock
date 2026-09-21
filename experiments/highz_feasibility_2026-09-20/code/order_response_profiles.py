#!/usr/bin/env python3
"""Profile the order offset without replacing nonlinear uncertainty by Hessian alone."""
import json
import numpy as np
from order_response_model import OUT,fit,SHA
summary=json.loads((OUT/'model_summary.json').read_text());grid=np.linspace(-.12,.12,13)
(OUT/'profile_protocol.json').write_text(json.dumps(dict(purpose='Conditional objective profiles after primary model fits, not additional hypothesis selection',grid_order_offset_km_s=grid.tolist(),models=['G1','A1'],transitions=[2600,2382],optimizer='Refit all remaining response, common velocities and row continua at each fixed offset; initialize from corresponding selected full fit',interpretation='No confidence calibration or Wilks claim; show actual bounded profile'),indent=2)+'\n')
results={}
for line in [2600,2382]:
 results[str(line)]={}
 for kind in ['G1','A1']:
  best=summary[str(line)]['full'][kind];rows=[]
  for i,delta in enumerate(grid):
   r=fit(f'profile_{line}_{kind}_{i:02d}',line,kind,start=best['parameters'],frozen={'order_offset':float(delta)})
   rows.append(dict(order_offset_m_s=1000*float(delta),chi2=r['chi2'],delta_chi2=r['chi2']-best['chi2'],optimizer_success=r['optimizer_success'],active_bounds=r['active_bounds'],source=r['name']+'.json'))
  results[str(line)][kind]=rows
  (OUT/'profile_summary.json').write_text(json.dumps(results,indent=2)+'\n')
(OUT/'profile_run_manifest.json').write_text(json.dumps(dict(code_sha256=SHA(__file__),model_summary_sha256=SHA(OUT/'model_summary.json'),protocol_sha256=SHA(OUT/'profile_protocol.json')),indent=2)+'\n')
