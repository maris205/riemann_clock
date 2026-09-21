#!/usr/bin/env python3
"""Isolated synthetic lines on real SQUAD grids with archived noise scales.

Observed flux is intentionally not fitted: this is a statistical feasibility
benchmark, not a line identification, physical fit or new-physics constraint.
"""
from pathlib import Path
import csv
import json
import hashlib
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
C_KMS = 299792.458
SEED = 20260921


def model(v, p):
    center, depth, width, offset, slope = p
    return offset + slope*v/25 - depth*np.exp(-.5*((v-center)/width)**2)


def jac(v, p):
    center, depth, width, _, _ = p
    dv = v-center
    exp = np.exp(-.5*(dv/width)**2)
    return np.column_stack([-depth*exp*dv/width**2, -exp,
                            -depth*exp*dv**2/width**3, np.ones_like(v), v/25])


def run():
    rng = np.random.default_rng(SEED)
    # Fixed illustrative post-instrument depth and Gaussian sigma. R is assumed,
    # NOT inferred from the actual coadd's changing line-spread function.
    resolving_power = 45000.
    width = np.sqrt((3./np.sqrt(2))**2+(C_KMS/resolving_power/2.354820045)**2)
    truth = np.array([.35, .25, width, 1., 0.])
    inputs = list(csv.DictReader((BASE/'data/processed/feii_coverage_audit.csv').open()))
    rows = []
    examples = []
    sources = {}
    for row in inputs:
        path = BASE/'data/processed'/f"{row['target']}_{row['transition']}_window.npz"
        sources[str(path.relative_to(BASE))] = hashlib.sha256(path.read_bytes()).hexdigest()
        data = np.load(path, allow_pickle=False)
        keep = data['valid'].astype(bool) & (np.abs(data['velocity_km_s']) < 25)
        v = data['velocity_km_s'][keep]
        # A window median avoids treating line-dependent photon errors as a
        # prediction for a different injected line depth.
        noise = float(np.median(data['error_normalized'][data['valid'].astype(bool)]))
        j = jac(v, truth)
        covariance = np.linalg.inv(j.T @ j/noise**2)
        centroid_sigma = float(np.sqrt(covariance[0,0]))
        adjacent = np.eye(len(v), k=1)+np.eye(len(v), k=-1)
        correlated_cov = noise**2*(np.eye(len(v))+.3*adjacent)
        correlated_sigma = float(np.sqrt(np.linalg.inv(j.T@np.linalg.solve(correlated_cov,j))[0,0]))
        snr = 1/noise
        eligible = snr >= 20 and len(v) >= 20
        out = {"target":row['target'], "transition":row['transition'],
               "z_abs_reference":float(row['z_abs_reference']), "N_grid_pixels_pm25kms":len(v),
               "median_continuum_SNR_assumed":snr,
               "Fisher_centroid_sigma_mps_diagonal":centroid_sigma*1000,
               "Fisher_centroid_sigma_mps_adjacent_rho03":correlated_sigma*1000,
               "MC_run":eligible, "n_MC":300 if eligible else 0,
               "MC_mean_bias_over_Fisher_sigma":None,
               "MC_std_over_Fisher_sigma":None,"MC_coverage95_fixed_Fisher":None,
               "MC_failed_or_bound_fits":None,
               "status":"SYNTHETIC_NOISE_BENCHMARK" if eligible else "LOW_SNR_LOCAL_FISHER_NOT_RELIABLE"}
        if eligible:
            centers=[]; bad=0
            for i in range(300):
                y = model(v,truth)+rng.normal(0,noise,len(v))
                fit=least_squares(lambda p:(model(v,p)-y)/noise,
                                  np.array([0., .2, width*1.1, 1., 0.]),
                                  jac=lambda p:jac(v,p)/noise,
                                  bounds=([-10,.001,.5,.5,-.5],[10,.95,12,1.5,.5]))
                bad+=int(not fit.success or np.any(fit.active_mask))
                centers.append(fit.x[0])
                if i==0 and row['transition'].startswith('FeII_1608'):
                    examples.append((row['target'],v,y,noise,fit.x))
            normalized=(np.array(centers)-truth[0])/centroid_sigma
            out.update(MC_mean_bias_over_Fisher_sigma=float(normalized.mean()),
                       MC_std_over_Fisher_sigma=float(normalized.std(ddof=1)),
                       MC_coverage95_fixed_Fisher=float(np.mean(np.abs(normalized)<norm.ppf(.975))),
                       MC_failed_or_bound_fits=bad)
        rows.append(out)
    result={"status":"SYNTHETIC_LINES; observed flux NOT fitted",
            "seed":SEED, "source_window_sha256":sources,
            "assumptions":{"independent_gaussian_noise_for_MC":True,"continuum_noise":"archived full-window median normalized error",
                           "real_grid":True,"assumed_R":resolving_power,"intrinsic_b_kms":3.,
                           "post_LSF_depth":.25,"post_LSF_sigma_kms":float(width),
                           "injected_center_kms":.35,"covariance_sensitivity_adjacent_rho":.3,
                           "nuisance_parameters_fitted":["depth","width","continuum_offset","continuum_slope"],
                           "excluded":["actual absorption profiles","blends","tellurics","real LSF","calibration bias","atomic wavelength uncertainty"]},
            "windows":rows}
    (BASE/'results/pixel_noise_pilot.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (BASE/'results/pixel_noise_pilot.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fig,axes=plt.subplots(1,2,figsize=(10,3.8),constrained_layout=True)
    for ax,(target,v,y,noise,p) in zip(axes,examples):
        ax.errorbar(v,y,yerr=noise,fmt='.',color='#467e97',alpha=.8,label='Synthetic pixels')
        vv=np.linspace(-25,25,300)
        ax.plot(vv,model(vv,truth),color='#252525',lw=1.5,label='Injected line')
        ax.plot(vv,model(vv,p),color='#b96344',ls='--',label='Recovered fit')
        ax.set(xlabel='Velocity (km/s)',ylabel='Normalized synthetic flux',title=f'{target}: Fe II 1608 grid')
        ax.legend(fontsize=8,frameon=False);ax.grid(alpha=.2)
    fig.suptitle('Mock lines / actual wavelength grids and median noise scales',fontsize=11)
    for suffix in ('png','pdf'):
        fig.savefig(BASE/f'figures/pixel_noise_pilot.{suffix}',dpi=180)
    plt.close(fig)
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    run()
