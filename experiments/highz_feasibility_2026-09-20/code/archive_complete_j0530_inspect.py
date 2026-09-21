"""Pre-fit source-profile inspection for the fixed J0530 campaign."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from archive_profile_pilot import atoms,C
OUT=ROOT/'results/archive_complete/J053007-250329';OUT.mkdir(parents=True,exist_ok=True)
d=np.load(ROOT/'data/processed/archive_expansion/J053007-250329.npz');a=atoms([1608,1611,2374,2382])
fig,axes=plt.subplots(4,2,figsize=(15,10),sharex='col',sharey=True)
for axrow,(key,aa) in zip(axes,a.items()):
 ref=np.sum(aa[:,0]*aa[:,1])/aa[:,1].sum();v=C*np.log(d['wavelength_AA']/(ref*3.141))
 for ax,lim in zip(axrow,[[-500,500],[-130,160]]):
  q=(v>=lim[0])&(v<=lim[1]);ax.plot(v[q],d['flux'][q],lw=.8);ax.fill_between(v[q],d['flux'][q]-d['expected_fluctuation'][q],d['flux'][q]+d['expected_fluctuation'][q],alpha=.16)
  bad=q&~d['valid'];ax.plot(v[bad],d['flux'][bad],'rx');ax.axhline(1,color='gray',lw=.5);ax.axhline(0,color='gray',lw=.5);ax.set_xlim(lim);ax.set_ylim(-.08,1.15);ax.grid(alpha=.2);ax.set_title(f'Fe II {key}');ax.axvline(-110,color='green',ls='--');ax.axvline(130,color='green',ls='--')
axes[-1,0].set_xlabel('Velocity relative to z = 2.141 (km/s)');axes[-1,1].set_xlabel('Velocity relative to z = 2.141 (km/s)')
fig.suptitle('J053007−250329: original coadd, before profile or relative-shift fitting\nGreen boundaries: candidate common window [−110, +130] km/s; bands: expected-fluctuation errors')
fig.tight_layout();fig.savefig(OUT/'pre_fit_profiles.png',dpi=160);fig.savefig(OUT/'pre_fit_profiles.pdf');plt.close(fig)
