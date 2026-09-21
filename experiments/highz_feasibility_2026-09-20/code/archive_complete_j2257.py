#!/usr/bin/env python3
"""Freeze and run the remaining J225719-100104 dominant-complex analysis."""
import argparse,json
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive_complete_model import ROOT,OUT,DATA,C,sha,run_campaign
from archive_profile_pilot import atoms
TARGET='J225719-100104';FOLDER=OUT/TARGET
def prepare():
 FOLDER.mkdir(parents=True,exist_ok=True)
 config={'target':TARGET,'z':1.836,'lines':[2344,2374,2382],'window':[-310.,110.],
  'component_bounds':[-290.,85.],'logtotal':15.0,'component_counts':[6,10,14,18,22],
  'centers':[-275,-241,-218,-184,-163,-151,-139,-130,-119,-108,-97,-85,-76,-65,-50,-40,-30,-18,20,33,55,63],
  'centers_by_n':{'6':[-241,-151,-108,-76,-30,60],
    '10':[-275,-241,-183,-150,-120,-95,-65,-30,22,60]},
  'seeds':[0,1,2],'max_nfev':350,'frozen_utc':datetime.now(timezone.utc).isoformat(),
  'selection_stage':'Frozen after raw profile inspection, before any conventional or relative-shift fit.',
  'selection_reason':'Dominant shared complex extends roughly -285 to +80km/s; common window provides margins. More blueward Fe2382 absorption around-330km/s is outside the modeled complex and has no obvious strong counterpart in this inspection; physical identity is not established.',
  'excluded_lines':{'1608':'Excluded by original catalogue-CNR policy; actualCNR~17 does not retroactively change primary line selection.',
    '2586,2600':'Retain original broad water-band/coverage exclusions.',
    '1611,2249,2260':'Not selected for primary multi-line fit; weak profiles are not a demonstrated clean unsaturated counterpart across this broad complex.'},
  'limitations':['Exploratory dominant-complex window, not proof the whole physical absorber is captured.',
    'Strong2382 saturation extends over150km/s; weaker2374 must constrain gas structure.',
    'Finite count/initialization search does not certify global optimum.'],
  'source_hashes':{'coadd':sha(DATA/f'{TARGET}.npz'),'metadata':sha(DATA/f'{TARGET}_metadata.json')},
 }
 path=FOLDER/'selection_frozen.json'
 if path.exists():raise FileExistsError('Refusing to overwrite previously frozen selection')
 path.write_text(json.dumps(config,indent=2)+'\n')
 d=np.load(DATA/f'{TARGET}.npz');ad=atoms([1608,2344,2374,2382,2586,2600]);fig,axes=plt.subplots(6,1,figsize=(13,10),sharex=True)
 for ax,(k,a) in zip(axes,ad.items()):
  ref=np.sum(a[:,0]*a[:,1])/a[:,1].sum();v=C*np.log(d['wavelength_AA']/(ref*(1+config['z'])));u=(abs(v)<550)&d['valid']
  ax.plot(v[u],d['flux'][u],lw=.6);ax.fill_between(v[u],d['flux'][u]-d['expected_fluctuation'][u],d['flux'][u]+d['expected_fluctuation'][u],alpha=.2)
  ax.axvspan(*config['window'],color='#d9ead3',alpha=.2);ax.set_ylim(-.1,1.2);ax.set_ylabel(str(k));ax.axhline(1,color='.4',lw=.6);ax.grid(alpha=.2)
 axes[-1].set_xlabel('Velocity relative to rounded catalogue redshift(km/s)');fig.suptitle(TARGET+' z1.836: pre-fit selection; shaded common window');fig.tight_layout();fig.savefig(FOLDER/'prefit_wide_inspection.png',dpi=140);fig.savefig(FOLDER/'prefit_wide_inspection.pdf');plt.close(fig)
 return path
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--prepare-only',action='store_true');args=a.parse_args()
 path=FOLDER/'selection_frozen.json'
 if not path.exists():prepare()
 if not args.prepare_only:print(json.dumps(run_campaign(path,workers=3),indent=2))
