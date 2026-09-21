#!/usr/bin/env python3
"""Freeze J0643 profile-selection metadata before fitting any relative shift."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from archive_profile_pilot import atoms,C
TARGET='J064326-504112'
D=ROOT/'data/processed/archive_expansion';OUT=ROOT/'results/archive_complete'/TARGET

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 s=np.load(D/f'{TARGET}.npz');meta=json.loads((D/f'{TARGET}_metadata.json').read_text())
 keys=[1608,1611,2249,2260,2344,2374,2382,2586,2600];atomic=atoms(keys)
 report={"target":TARGET,"z":2.659,"common_window_km_s":[-145,100],"recommendation_primary_lines":[1608,1611,2374,2382],"decision_basis":"Raw profiles inspected at ±500 km/s, prior to fitting shifts. Shared main complex from approximately −135 to +60 km/s; weak 1611 constrains saturated gas. Main complex isolated with fixed common window. No residual clipping.","separate_features":"A feature near −157 km/s in 2382 lies outside the main window; corresponding much weaker 1608 structure may exist, so its blend identity is not established. 1611 has separate absorption near −170 and +133 km/s outside the window. 1608 has strong extra absorption near +340 km/s. These are not treated as evidence for anomalous relative shifts.","spectra_sha256":sha(D/f'{TARGET}.npz'),"metadata_sha256":sha(D/f'{TARGET}_metadata.json'),"atomic_sha256":sha(ROOT/'data/raw/espresso_null/MM_VPFIT_2013-11-10.dat'),"boundary_limitation":"Window left boundary lies between the nearby blue feature and main complex, with only short continuum separation. An unmodelled wing could affect precision; no blind or anomaly-independent global sample claim.","selection_conditioning":"Window and line choices use observed morphology and are exploratory. They precede all relative-shift fits to this target.","lines":[]}
 fig,axes=plt.subplots(9,1,figsize=(15,20),sharex=True)
 for ax,k in zip(axes,keys):
  a=atomic[k];ref=float(np.sum(a[:,0]*a[:,1])/sum(a[:,1]));v=C*np.log(s['wavelength_AA']/(ref*3.659));show=(abs(v)<500)&s['valid'];u=(v>=-145)&(v<=100);good=u&s['valid']&np.isfinite(s['expected_fluctuation'])&(s['expected_fluctuation']>0)&np.isfinite(s['flux'])
  ax.step(v[show],s['flux'][show],where='mid',lw=.8);ax.plot(v[show],s['expected_fluctuation'][show],c='orange',alpha=.6);ax.axhline(1,c='gray',lw=.5);ax.axvline(-145,c='red',lw=.8);ax.axvline(100,c='red',lw=.8);ax.set_ylim(-.15,1.5);ax.text(.01,.12,str(k),transform=ax.transAxes)
  report['lines'].append(dict(line=k,reference_AA=ref,ndata=int(good.sum()),window_total=int(u.sum()),observed_wavelength_range_AA=[float(s['wavelength_AA'][u].min()),float(s['wavelength_AA'][u].max())],median_expected_error=float(np.median(s['expected_fluctuation'][good])),median_stat_error=float(np.median(s['error'][good])),min_flux=float(np.min(s['flux'][good])),below_0p1_pixels=int(np.sum(s['flux'][good]<.1)),expected_valid_source_indices=np.flatnonzero(good).tolist()))
 axes[-1].set_xlabel('Velocity (km/s); isotope-weighted laboratory reference at z=2.659');fig.tight_layout();fig.savefig(OUT/'selection_full_context.png',dpi=130);fig.savefig(OUT/'selection_full_context.pdf');plt.close(fig)
 (OUT/'pre_fit_selection.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
