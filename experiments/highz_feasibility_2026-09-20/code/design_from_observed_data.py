#!/usr/bin/env python3
"""Data-informed selection and design arithmetic; no cosmic parameter fit."""
from pathlib import Path
import csv
import json
import numpy as np

BASE=Path(__file__).resolve().parents[1]


def run():
    quality=list(csv.DictReader((BASE/'data/processed/feii_coverage_audit.csv').open()))
    forecast=json.loads((BASE/'results/forecast_summary.json').read_text())
    times={r['z']:r for r in forecast['time_leverage']}
    g_low=times[1.1508]['G2_z3_normalized'];g_high=times[3.025]['G2_z3_normalized']
    rows=[]
    for r in quality:
        snr=float(r['median_continuum_to_error_per_native_pixel'])
        line=int(float(r['rest_wavelength_AA']))
        high=r['target']=='J034943-381030'
        flag=('LOW_SNR_EXCLUDE' if high and line in (2586,2600) else
              'TELLURIC_SCREEN_ONLY' if high and line in (2344,2374,2382) else
              'PUBLISHED_BLEND_REVIEW' if not high and line in (2344,2586) else
              'DEEP_SATURATION_REVIEW' if not high and line in (2382,2600) else
              'CANDIDATE_REQUIRES_FULL_CALIBRATION')
        rows.append({'target':r['target'],'transition':r['transition'],
             'z_abs':float(r['z_abs_reference']),'observed_wavelength_A':float(r['nominal_observed_wavelength_AA']),
             'median_CNR_per_1p3kms_pixel':snr,'design_status':flag,
             'photon_only_exposure_factor_to_CNR30':max(1.,(30/snr)**2),
             'photon_only_exposure_factor_to_CNR50':max(1.,(50/snr)**2),
             'exposure_factor_caveat':'Assumes identical throughput, background regime and independent photon noise; not telescope time or correction for blends/tellurics'})
    # Even an optimistic six-line low-z measurement + one high-z transition
    # cannot separate B from the response-mode intercept and absorber redshifts.
    # Toy K is only a numerical rank demonstration. No claim of atomic K.
    k=np.array(forecast['synthetic_response']['K'])
    wave=np.array(forecast['synthetic_response']['rest_vacuum_wavelength_grid_A_lookup_only'])
    x=np.log(wave);x=(x-x.mean())/np.ptp(x)
    va_low=np.r_[np.ones(6),0.];va_high=np.r_[np.zeros(6),1.]
    slope=np.r_[x,0.];kall=np.r_[k,k[0]]
    signal=np.r_[g_low*k,g_high*k[0]]
    nuisance=np.column_stack([va_low,va_high,slope,kall])
    full=np.column_stack([nuisance,signal])
    residual=signal-nuisance@np.linalg.lstsq(nuisance,signal,rcond=None)[0]
    null=np.array([0.,-(g_high-g_low)*k[0],0.,-g_low,1.])
    assert np.linalg.norm(full@null)<1e-12
    assert np.linalg.norm(residual)<1e-12
    summary={'status':'DESIGN_DERIVED_FROM_OBSERVED_COVERAGE; no cosmic estimate',
        'quality_rows':rows,
        'redshift_leverage':{'z_low':1.1508,'z_high':3.025,'G2_low':g_low,'G2_high':g_high,
                            'delta_G2':g_high-g_low,
                            'B100mps_mode_difference_mps':100*(g_high-g_low),
                            'B30mps_mode_difference_mps':30*(g_high-g_low),
                            'note':'B is an unspecified hypothetical response amplitude at z=3, not measured or predicted'},
        'optimistic_identifiability':{'rows':7,'low_z_lines_assumed':6,'high_z_lines_assumed':1,
                            'nuisance_rank':int(np.linalg.matrix_rank(nuisance)),
                            'full_rank':int(np.linalg.matrix_rank(full)),
                            'full_columns':5,'full_singular_values':np.linalg.svd(full,compute_uv=False).tolist(),
                            'projected_signal_norm':float(np.linalg.norm(residual)),
                            'analytic_null_vector':null.tolist(),
                            'reason':'One high-z line is absorbed by its free redshift; at one remaining epoch B and b0 are proportional. Independent arbitrary calibration priors cannot be assumed.',
                            'K_is_artificial':True},
        'prior_forecast_status':'24/60-sightline six-clean-line forecasts remain prospective, not achieved by these two spectra'}
    with (BASE/'results/data_informed_design.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (BASE/'results/data_informed_design.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'redshift_leverage':summary['redshift_leverage'],
                      'identifiability':summary['optimistic_identifiability'],
                      'low_SNR_exposure_factors':[r for r in rows if r['design_status']=='LOW_SNR_EXCLUDE']},indent=2))


if __name__=='__main__':run()
