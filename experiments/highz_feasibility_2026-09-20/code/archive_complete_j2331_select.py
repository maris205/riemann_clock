#!/usr/bin/env python3
"""Freeze J233156-090802 profile selection before any relative-shift fit."""
from pathlib import Path
from datetime import datetime, timezone
import json, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive_profile_pilot import atoms, C
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/archive_complete/J233156-090802'
DATA=ROOT/'data/processed/archive_expansion/J233156-090802.npz'
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 a=np.load(DATA);keys=[1608,1611,2374,2382];at=atoms(keys)
 fig,axes=plt.subplots(4,1,figsize=(15,11),sharex=True)
 pixels={}
 for key,ax in zip(keys,axes):
  atomic=at[key];ref=np.sum(atomic[:,0]*atomic[:,1])/atomic[:,1].sum()
  v=C*np.log(a['wavelength_AA']/(ref*3.143))
  show=(abs(v)<500)&a['valid'];ax.plot(v[show],a['flux'][show],lw=.7)
  ax.fill_between(v[show],a['flux'][show]-a['expected_fluctuation'][show],a['flux'][show]+a['expected_fluctuation'][show],alpha=.2)
  ax.axhline(1,color='gray',lw=.5);ax.axvspan(-480,500,color='green',alpha=.05)
  ax.axvline(-480,color='green',lw=.6);ax.axvline(500,color='green',lw=.6)
  ax.set_ylabel(f'Fe II {key}');ax.grid();ax.set_ylim(-.1,1.25)
  within=(v>=-480)&(v<=500)
  use=within&a['valid']&np.isfinite(a['flux'])&np.isfinite(a['expected_fluctuation'])&(a['expected_fluctuation']>0)
  pixels[str(key)]={'within_window':int(within.sum()),'source_valid_expected_noise':int(use.sum()),'native_dv_km_s':float(np.median(np.diff(v)))}
 axes[-1].set_xlabel('Velocity relative to rounded catalogue z = 2.143 (km/s)')
 fig.suptitle('J233156−090802: expanded-complex selection inspected before shift fitting; green interval [−480, +500] km/s')
 fig.tight_layout();fig.savefig(OUT/'full_500kms_inspection.png',dpi=150);fig.savefig(OUT/'full_500kms_inspection.pdf');plt.close(fig)
 policy={'target':'J233156-090802','z':2.143,'lines':[1608,2374,2382],'window':[-480.,500.],'component_bounds':[-465.,480.],
         'centers':[-425.,-400.,-375.,-355.,-335.,-315.,-295.,-270.,-240.,-220.,-200.,-150.,-115.,-80.,25.,72.,180.,360.,455.],
         'logtotal':15.4,'component_counts':[6,10,14,18,22,26],'seeds':[0,1,2],'max_nfev':350,
         'selection_stage':'Frozen after raw profile and AOD inspection, before nonlinear or relative-shift fitting.',
         'selection_reason':'Shared Fe absorption occupies approximately -430 to +85 km/s. The common window covers the dominant visible shared complex with continuum margins; detached red features are not assumed absent in Fe. The saturated strong line must be constrained jointly by 1608 and 2374.',
         'excluded_lines':{'1611':'Weak absorption cannot constrain the complex gas independently and shows unrelated blue absorption; not in original three-line common-detection set.','2249,2260,2344':'Retain earlier broad atmospheric-water warnings.','2586,2600':'Retain earlier coverage/quality exclusions.'},
         'outside_window_features':'Detached troughs strongest in 2382 near +180,+360,+450 km/s lie outside the chosen dominant-complex window; weak Fe counterparts and external blends are both possible. Their identification is not established. No internal pixels clipped by fitted residuals.',
         'limitations':['Exploratory data-inspected selection, not blind confirmation.','Strong 2382 saturated across roughly 400 km/s.','Twenty-six shared gas components may still be insufficient; inadequacy is retained as an outcome.','Apparent 1608 continuum near 0.92 requires continuum flexibility.','No alpha or cosmic-time-law parameter is fitted.'],
         'pixel_diagnostics':pixels,'source_sha256':sha(DATA),'selection_script_sha256':sha(__file__),'frozen_utc':datetime.now(timezone.utc).isoformat()}
 p=OUT/'selection.json'
 if p.exists():
  old=json.loads(p.read_text());assert all(old[k]==policy[k] for k in ['target','z','lines','window','component_bounds','centers','component_counts','seeds','source_sha256']), 'Frozen scientific selection changed.'
 else:p.write_text(json.dumps(policy,indent=2)+'\n')
 print(json.dumps(pixels,indent=2))
if __name__=='__main__':main()
