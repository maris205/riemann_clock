#!/usr/bin/env python3
"""Actual native-pixel FeII quality and apparent-optical-depth archive audit.

No line displacement is fitted and no cosmic-time function is fitted. Fixed
catalogue-redshift windows are readiness diagnostics, not precision velocities.
"""
from __future__ import annotations
import csv, hashlib, json, re
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from archive_screen import ROOT, OUT, LINES, C
from archive_strict_screen import BANDS

DATA=ROOT/'data/processed/archive_expansion'
RAW=ROOT/'data/raw/archive_expansion'

QUALITY_POLICY={
 'full_half_width_km_s':250.,'core_half_width_km_s':100.,
 'minimum_valid_fraction_full':.98,'minimum_native_CNR':15.,
 'minimum_core_EW_diagonal_SNR':5.,'minimum_detected_strong_lines':3,
 'weak_line_requirement':['2374','2586'],
 'common_grid_spacing_km_s':5.,'common_absorption_diagonal_SNR':3.,
 'minimum_common_grid_points_at_least_3_strong_lines':3,
 'weak_line_minimum_unsaturated_detection_pixels':5,
 'common_grid_note':'Interpolated per-pixel nominal absorption S/N on a 5-km/s grid; not independent bin-averaged significance. Both bracketing native pixels must be valid.',
 'AOD_flux_floor':'max(0.02,3*normalized_error); floor-affected pixels are lower-limit diagnostics, not recovered intrinsic optical depths',
 'AOD_scaling':'Na(v)=3.768e14*tau_a/(f*lambda_AA) cm^-2/(km/s)',
 'eligibility_interpretation':'Coverage, S/N and multiplet-detection readiness only; no atmospheric template, sky, blending, wavelength-calibration, covariance or physical-gas-model certification.',
}

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def scalar(x): return float(x) if np.isfinite(x) else None
def overlap(lo,hi,bands):return any(lo<b and hi>a for a,b in bands)
def longest_false(a):
    bad=np.r_[False,~a,False];edges=np.flatnonzero(bad[1:]!=bad[:-1])
    return int(max(edges[1::2]-edges[::2],default=0))

def analyze_line(target,z,l,arrays,meta):
    wave=arrays['wavelength_AA'];v=C*np.log(wave/(l['rest_AA']*(1+z)))
    ii=np.flatnonzero(abs(v)<=250);vv=v[ii];f=arrays['flux'][ii];e=arrays['error'][ii];g=arrays['valid'][ii]
    dv=C*np.log(10)*meta['header']['CD1_1'];core=(abs(vv)<=100)&g;side=(abs(vv)>=150)&g
    # Use actual log-grid pixel-edge widths, not a constant wavelength approximation.
    dl=wave[ii]*(np.exp(dv/2/C)-np.exp(-dv/2/C))/(1+z)
    lo,hi=np.array(l['window_AA'])*np.exp(np.array([-30.,30.])/C)
    water=overlap(lo,hi,BANDS)
    sky=overlap(lo,hi,[(x-1,x+1) for x in [5578.,5891.,5897.]])
    artifact=overlap(lo,hi,[(4713,4719),(4741,4747),(5237,5243),(5580,5800)])
    ew=float(np.sum((1-f[core])*dl[core]));es=float(np.sqrt(np.sum((e[core]*dl[core])**2)))
    ewfull=float(np.sum((1-f[g])*dl[g]));esfull=float(np.sqrt(np.sum((e[g]*dl[g])**2)))
    cnr=float(np.median(1/e[g])) if g.any() else np.nan
    ng=int(g.sum());n=len(ii);frac=ng/n if n else 0
    floor=np.maximum(.02,3*np.where(e>0,e,np.inf));tau=np.full(n,np.nan)
    tau[g]=-np.log(np.maximum(f[g],floor[g]));na=tau*3.768e14/(l['f']*l['rest_AA'])
    sat=core&(f<=np.maximum(.1,3*e))
    useful=core&(f>np.maximum(.1,3*e))&(f<.95)&((1-f)>3*e)
    reasons=[]
    if not l['eligible']:reasons.extend(l['rejections'])
    if frac<.98:reasons.append('native_valid_fraction_below_0p98')
    if not np.isfinite(cnr) or cnr<15:reasons.append('native_CNR_below_15')
    if es<=0 or ew/es<5:reasons.append('core_EW_diagonal_detection_proxy_below_5')
    # Fixed edge diagnostic flags a potentially truncated absorption complex.
    edge=g&(abs(vv)>200)
    edge_sig=np.sum((1-f[edge])*dl[edge])/np.sqrt(np.sum((e[edge]*dl[edge])**2)) if edge.any() else np.nan
    near=min(meta['tables']['1']['rows'],key=lambda r:abs(r['Wavelength']-l['observed_AA']))
    row={'target':target,'z_abs_catalogue':z,'line':l['line'],'rest_AA':l['rest_AA'],'f':l['f'],
      'observed_AA_catalogue':l['observed_AA'],'metadata_eligible':l['eligible'],
      'native_pixels':n,'valid_pixels':ng,'valid_fraction':frac,'longest_invalid_run_pixels':longest_false(g),
      'actual_log_grid_spacing_km_s':float(dv),
      'native_CNR':scalar(cnr),'sideband_median_flux':scalar(np.median(f[side])) if side.any() else None,
      'core_EW_rest_AA':ew,'core_EW_diagonal_error_AA':es,'core_EW_diagonal_SNR':ew/es if es else None,
      'full_EW_rest_AA':ewfull,'full_EW_diagonal_error_AA':esfull,
      'core_to_full_EW_ratio':ew/ewfull if ewfull>0 else None,
      'core_low_flux_pixels':int(sat.sum()),'core_valid_pixels':int(core.sum()),
      'core_low_flux_fraction':float(sat.sum()/core.sum()) if core.any() else None,
      'core_resolved_unsaturated_detection_pixels':int(useful.sum()),
      'core_AOD_integral_cm2':float(np.nansum(na[core])*dv),
      'full_window_edge_EW_diagonal_SNR':scalar(edge_sig),
      'broader_water_band_warning':water,'strong_sky_feature_warning':sky,'reported_instrument_artifact_warning':artifact,
      'nearest_nominal_R':float(near['NomResolPower']),'nearest_arc_R':float(near['ArcResolPower']),
      'median_contributing_extracted_pixels':float(np.median(arrays['n_contrib_after'][ii][g])) if g.any() else None,
      'detected_primary_candidate':not reasons,'detected_strict_candidate':not reasons and not water,
      'rejections':reasons}
    ident=f'{target}_z{z:.3f}_FeII{l["line"]}'
    np.savez_compressed(DATA/(ident+'_window.npz'),source_indices=ii,velocity_km_s=vv,wavelength_AA=wave[ii],
       flux=f,error=e,valid=g,AOD_tau=tau,AOD_Na_per_km_s=na,core_mask=core,
       uncertainty_floor_affected=g&(f<=floor))
    return row,(vv,f,e,g,na)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    policy_file=OUT/'actual_quality_policy.json'
    if not policy_file.exists():policy_file.write_text(json.dumps({'frozen_utc':datetime.now(timezone.utc).isoformat(),**QUALITY_POLICY},indent=2)+'\n')
    screen=json.loads((OUT/'catalogue_screen.json').read_text())
    rows=[];absorbers=[];figcache={}
    with PdfPages(OUT/'all_absorber_profiles.pdf') as pdf:
      for r in screen['all_absorbers']:
        if not r['metadata_eligible']:continue
        target=r['target'];z=r['z_abs'];path=DATA/(target+'.npz')
        if not path.exists():
            absorbers.append({'target':target,'z_abs':z,'status':'download_missing'});continue
        arrays=np.load(path);meta=json.loads((DATA/(target+'_metadata.json')).read_text())
        line_rows=[];profiles={}
        for l in r['lines']:
            row,p=analyze_line(target,z,l,arrays,meta);rows.append(row);line_rows.append(row);profiles[l['line']]=p
        primary=[l for l in line_rows if l['f']>=.01 and l['detected_primary_candidate']]
        strict=[l for l in line_rows if l['f']>=.01 and l['detected_strict_candidate']]
        weak=[l['line'] for l in strict if l['line'] in ['2374','2586'] and l['core_resolved_unsaturated_detection_pixels']>=5]
        grid=np.arange(-100,100.1,5);detection=[]
        for l in strict:
            v,f,e,g,na=profiles[l['line']]
            ix=np.searchsorted(v,grid);inside=(ix>0)&(ix<len(v));ix=np.clip(ix,1,len(v)-1)
            bracket_valid=inside&g[ix-1]&g[ix]
            detection.append(bracket_valid&(np.interp(grid,v[g],(1-f[g])/e[g],left=0,right=0)>=3))
        common=int(np.sum(np.sum(detection,axis=0)>=3)) if detection else 0
        possible=len(strict)>=3 and bool(weak) and common>=3
        upl=(RAW/(target+'.upl')).read_text() if (RAW/(target+'.upl')).exists() else ''
        # Count recorded correction lines; do not infer achieved calibration accuracy.
        vsht=[line.strip() for line in upl.splitlines() if 'VSHT' in line]
        absrow={'target':target,'z_abs':z,'z_em':r['z_em'],'status':'actual_arrays_audited',
          'n_primary_detected_strong_lines':len(primary),'primary_detected_strong_lines':[l['line'] for l in primary],
          'n_strict_detected_strong_lines':len(strict),'strict_detected_strong_lines':[l['line'] for l in strict],
          'resolved_weak_line_candidates':weak,'common_5km_s_grid_points_with_3line_absorption':common,
          'passes_computational_readiness':possible,
          'readiness_label':'multiplet detected; model development candidate, calibration/blends/atmosphere unvalidated' if possible else 'fails fixed conservative readiness gate',
          'n_exposures':meta['n_exposures'],'native_dispersion_km_s':meta['native_dispersion_km_s'],
          'observation_dates':meta['observation_date_range'],'UPL_VSHT_records':len(vsht),
          'UPL_VSHT_records_excerpt':vsht[:3],
          'artifact_warning_candidate_lines':[l['line'] for l in strict if l['reported_instrument_artifact_warning']],
          'edge_absorption_warning_lines':[l['line'] for l in strict if (l['full_window_edge_EW_diagonal_SNR'] or 0)>5],
          'line_rows':line_rows}
        absorbers.append(absrow)
        fig,axes=plt.subplots(3,3,figsize=(11.7,8.3),sharex=True,sharey=True)
        for ax,l in zip(axes.flat,line_rows):
            v,f,e,g,na=profiles[l['line']]
            ax.plot(v[g],f[g],lw=.6,color='#253b52');ax.fill_between(v[g],f[g]-e[g],f[g]+e[g],color='#7fa5c1',alpha=.25,lw=0)
            ax.axhline(1,color='.55',lw=.5);ax.axvline(-100,color='.75',ls=':',lw=.6);ax.axvline(100,color='.75',ls=':',lw=.6)
            label='candidate' if l['detected_strict_candidate'] else ('H2O warning' if l['detected_primary_candidate'] and l['broader_water_band_warning'] else 'excluded / weak')
            ax.set_title(f'Fe II {l["line"]}: {label}\nCNR {l["native_CNR"]:.0f}' if l['native_CNR'] is not None else f'Fe II {l["line"]}: no valid pixels',fontsize=9)
            ax.set_xlim(-250,250);ax.set_ylim(-.12,1.22);ax.grid(alpha=.15)
        fig.suptitle(f'{target}, catalogue z_abs={z:.3f}; {len(strict)} detected lines outside broad water bands\nNative coadd pixels; catalogue zero is rounded. No relative displacement or cosmic-time fit.',fontsize=12)
        fig.supxlabel('Velocity relative to rounded catalogue redshift (km/s)');fig.supylabel('Continuum-normalized flux')
        fig.tight_layout(rect=[.02,.02,1,.92]);pdf.savefig(fig)
        if possible or target in [s['target'] for s in screen['selected']]:fig.savefig(OUT/f'{target}_z{z:.3f}_profiles.png',dpi=130)
        plt.close(fig)
        if possible:figcache[(target,z)]=(line_rows,profiles)
        print(target,z,'strict',len(strict),'weak',weak,'common',common,'ready',possible,flush=True)
    with (OUT/'window_quality.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (OUT/'absorber_quality.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),
       'policy_sha256':digest(policy_file),'results':absorbers},indent=2)+'\n')
    # AOD overlays expose saturation/shape disagreement; these are not velocity fits.
    with PdfPages(OUT/'candidate_AOD_diagnostics.pdf') as pdf:
      for (target,z),(line_rows,profiles) in figcache.items():
        fig,ax=plt.subplots(figsize=(10,5))
        for l in line_rows:
          if l['f']>=.01 and l['detected_strict_candidate']:
            v,f,e,g,na=profiles[l['line']];ax.plot(v[g],na[g]/1e12,label=f'Fe II {l["line"]}',lw=.8)
        ax.set(xlim=(-150,150),ylim=(-.5,None),xlabel='Velocity relative to rounded catalogue redshift (km/s)',ylabel=r'Apparent column density ($10^{12}$ cm$^{-2}$ / km s$^{-1}$)',
          title=f'{target}, z_abs={z:.3f}: AOD consistency diagnostic\nLow flux uses a 3-error floor; saturation, blends and LSF remain unresolved')
        ax.legend();ax.grid(alpha=.2);fig.tight_layout();pdf.savefig(fig)
        fig.savefig(OUT/f'{target}_z{z:.3f}_AOD.png',dpi=140);plt.close(fig)
    summary={'created_utc':datetime.now(timezone.utc).isoformat(),
       'candidate_absorbers':len(absorbers),'downloaded_unique_sightlines':len({a['target'] for a in absorbers if a['status']=='actual_arrays_audited'}),
       'actual_windows':len(rows),'computational_readiness_count':sum(a.get('passes_computational_readiness',False) for a in absorbers),
       'readiness_candidates':[{k:a[k] for k in ['target','z_abs','strict_detected_strong_lines','resolved_weak_line_candidates','common_5km_s_grid_points_with_3line_absorption','edge_absorption_warning_lines','artifact_warning_candidate_lines']} for a in absorbers if a.get('passes_computational_readiness')],
       'no_inference':['No physical line shifts estimated','No alpha measurement','No time-law fit','No calibrated precision-constant measurement certified'],
       'hashes':{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),ROOT/'data/raw/espresso_null/MM_VPFIT_2013-11-10_noiso.dat',policy_file,OUT/'absorber_quality.json',OUT/'window_quality.csv']}}
    (OUT/'quality_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
