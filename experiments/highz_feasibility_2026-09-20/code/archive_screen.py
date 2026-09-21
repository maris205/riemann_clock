#!/usr/bin/env python3
"""Fixed, metadata-only screen of the published SQUAD DLA catalogue.

No line offsets enter the selection. Ranking is a convenience allocation of
archive processing, not a representative population or cosmology sample.
"""
from __future__ import annotations
import csv, hashlib, json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/archive_expansion'
C = 299792.458
LINES = {
    '1608': (1608.450852, .0577), '1611': (1611.200369, .00138),
    '2249': (2249.875472, .00182), '2260': (2260.779108, .00244),
    '2344': (2344.212747, .114), '2374': (2374.460064, .0313),
    '2382': (2382.763995, .32), '2586': (2586.649312, .0691),
    '2600': (2600.172114, .239),
}
POLICY = {
    'version': 1,
    'frozen_before_downloading_additional_spectra': True,
    'half_window_km_s': 250.,
    'metadata_CNR_per_2p5_km_s_min': 20.,
    'minimum_strong_lines_f_ge_0p01': 3,
    'quasar_proximity_km_s': 3000.,
    'red_of_quasar_Lyalpha_km_s': 3000.,
    'observed_wavelength_range_AA': [4000., 9000.],
    'conservative_telluric_exclusion_AA': [[6270., 6330.], [6860., 6960.],
       [7160., 7340.], [7580., 7700.], [8100., 8400.], [8900., 10000.]],
    'telluric_note': 'Conservative analyst-defined bands, not a measured transmission mask. Passing does not prove freedom from weak tellurics, sky or unrelated absorption.',
    'previous_targets_excluded': ['J051707-441055', 'J034943-381030'],
    'rank': 'third-highest proxy CNR among eligible strong lines, then count; descending. Top 5 distinct sightlines; retain all rejection reasons.',
    'actual_pixel_screen': {'minimum_valid_fraction': .98,
       'minimum_median_inverse_error_per_native_pixel': 15.,
       'minimum_absorption_equivalent_width_significance_diagonal': 5.,
       'minimum_strong_detected_lines': 3},
    'profile_diagnostics': 'Full fixed +/-250 km/s windows, then predeclared +/-100 km/s core diagnostic. No relative-velocity based rejection; no sub-km/s or cosmic-time inference from moment centroids.',
}

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    quasars = list(csv.DictReader((ROOT/'data/raw/DR1_quasars_master.csv').open()))
    lookup = {r['Name_Adopt']: r for r in quasars}
    dlas = list(csv.DictReader((ROOT/'data/raw/DR1_DLAs.csv').open()))
    rows = []
    for r in dlas:
        q = lookup[r['Name_Adopt']]
        cover = [tuple(map(float, s.split('-'))) for s in q['WavCoverage'].split(',') if s]
        cnr = list(map(float, q['CNR'].split(',')))
        for zstr in r['DLAzabs'].split(','):
            z = float(zstr); zem = float(q['zem_Adopt']); lineinfo = []
            for key,(rest,f) in LINES.items():
                wl = rest*(1+z)
                lo,hi = wl*np.exp(np.array([-1,1])*POLICY['half_window_km_s']/C)
                proxy = cnr[int(np.argmin(np.abs(np.array([3500,4500,5500,6500,7500])-wl)))]
                reasons = []
                if not any(a<=lo and hi<=b for a,b in cover): reasons.append('catalogue_coverage_gap')
                if lo<4000 or hi>9000: reasons.append('outside_fixed_wavelength_range')
                if any(lo<b and hi>a for a,b in POLICY['conservative_telluric_exclusion_AA']): reasons.append('conservative_telluric_band')
                if lo <= 1215.6701*(1+zem)*np.exp(3000/C): reasons.append('forest_or_Lyalpha_proximity')
                if proxy<20: reasons.append('low_catalogue_CNR_proxy')
                lineinfo.append({'line':key,'rest_AA':rest,'f':f,'observed_AA':wl,'CNR_proxy':proxy,
                    'window_AA':[lo,hi],'eligible':not reasons,'rejections':reasons})
            eligible = [x for x in lineinfo if x['eligible'] and x['f']>=.01]
            reasons = []
            if len(eligible)<3: reasons.append('fewer_than_3_eligible_strong_lines')
            if C*np.log((1+zem)/(1+z))<3000: reasons.append('absorber_quasar_proximity')
            if q['Name_Adopt'] in POLICY['previous_targets_excluded']: reasons.append('already_analyzed_target')
            rows.append({'target':q['Name_Adopt'],'z_abs':z,'z_em':zem,'DLA_reference':r['DLAReference'],
                'metadata_eligible':not reasons,'rejections':reasons,'lines':lineinfo,
                'n_strong_eligible':len(eligible),'rank_score':sorted([x['CNR_proxy'] for x in eligible],reverse=True)[2] if len(eligible)>=3 else 0,
                'exposure_s':float(q['ExpTime'] or 0),'dispersion_km_s':float(q['Dispersion'] or 0),
                'catalogue_Spec_status':q['Spec_status']})
    ranked = sorted([r for r in rows if r['metadata_eligible']],key=lambda r:(-r['rank_score'],-r['n_strong_eligible'],r['target'],r['z_abs']))
    chosen=[]; seen=set()
    for r in ranked:
        if r['target'] not in seen:
            chosen.append({'target':r['target'],'z_abs':r['z_abs'],'rank_score':r['rank_score'],'n_strong_eligible':r['n_strong_eligible']}); seen.add(r['target'])
        if len(chosen)==5: break
    result = {'created_utc':datetime.now(timezone.utc).isoformat(),'policy':POLICY,
        'source_hashes':{str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'data/raw/DR1_quasars_master.csv',ROOT/'data/raw/DR1_DLAs.csv',Path(__file__)]},
        'master_rows':len(quasars),'final_spectra_positive_exposure_count':sum(float(r['NumExp'] or 0)>0 for r in quasars),
        'DLA_sightlines':len(dlas),'DLA_absorbers':len(rows),'metadata_eligible_absorbers':len(ranked),
        'selected':chosen,'ranked_absorbers':[{'target':r['target'],'z_abs':r['z_abs'],'rank_score':r['rank_score'],'n_strong_eligible':r['n_strong_eligible']} for r in ranked],
        'all_absorbers':rows}
    (OUT/'catalogue_screen.json').write_text(json.dumps(result,indent=2)+'\n')
    (OUT/'selection_policy.json').write_text(json.dumps({'frozen_utc':result['created_utc'],**POLICY},indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['master_rows','final_spectra_positive_exposure_count','DLA_sightlines','DLA_absorbers','metadata_eligible_absorbers','selected']},indent=2))
if __name__=='__main__': main()
