"""Prospective conditional predictions; no measured spectrum or offset is fitted."""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT.parent / 'highz_feasibility_2026-09-20'
H0, OM, TSTAR = 67.4, .315, 5.391247e-44
MPC_KM, YEAR = 3.0856775814913673e19, 31557600.
Z = np.array([1., 1.5, 2.])


def age_gyr(z):
    z = np.asarray(z)
    return 2 * MPC_KM / H0 / (3 * np.sqrt(1-OM)) * np.arcsinh(
        np.sqrt((1-OM)/OM) / (1+z)**1.5) / (YEAR*1e9)


def f(z, p=2):
    l0 = np.log(age_gyr(0)*YEAR*1e9/TSTAR)
    return np.expm1(-p*np.log1p(np.log(age_gyr(z)/age_gyr(0))/l0))


def h(z, p=2):
    """Endpoint-normalized time law H(z_low)=0, H(z_high)=1."""
    return (f(z,p)-f(Z[0],p))/(f(Z[-1],p)-f(Z[0],p))


def power(effect, sigma, alpha=.01):
    critical = norm.isf(alpha/2)
    return norm.sf(critical-effect/sigma) + norm.cdf(-critical-effect/sigma)


def generate():
    for directory in ['results','figures']:
        (ROOT/directory).mkdir(exist_ok=True)
    ratio = float(h(Z[1]))
    contrast = np.array([-(1-ratio),1.,-ratio])
    # The 100 m/s endpoint contrast is an arbitrary illustration, NOT a fit.
    example_B = 100.
    alternatives = {}
    for p in [1,3]:
        alternatives[f'p{p}'] = dict(
            middle_fraction=float(h(Z[1],p)),
            middle_difference_at_example_B_m_s=float(example_B*(h(Z[1],p)-ratio)))
    scenarios = []
    threshold = brentq(lambda snr: power(snr,1.)-.9, 0.,10.)
    for sigma_single,tau_group in [(50.,10.),(100.,20.),(200.,30.)]:
        denominator = (example_B/threshold)**2/2-tau_group**2
        n90 = int(np.ceil(sigma_single**2/denominator)) if denominator>0 else None
        points = []
        for n in [4,8,12,24,48]:
            se = float(np.sqrt(2*(sigma_single**2/n+tau_group**2)))
            points.append(dict(n_per_endpoint_group=n, endpoint_difference_sigma_m_s=se,
                               power_if_B_100=float(power(example_B,se)),
                               effect_for_90pct_power_m_s=float(threshold*se)))
        scenarios.append(dict(single_target_random_sigma_m_s=sigma_single,
                              nonaveraging_group_sigma_m_s=tau_group,
                              required_n_per_endpoint_for_90pct_power_at_B100=n90,
                              maximum_power_at_infinite_n=float(power(example_B,np.sqrt(2)*tau_group)),
                              examples=points))
    atomic_path=OLD/'data/atomic/selected_transitions.csv'
    atomic=list(csv.DictReader(atomic_path.open()))
    wavelengths={row['label']:float(row['ritz_vacuum_wavelength_A']) for row in atomic
                 if row['label'] in ['Fe II 2260','Fe II 2344','Fe II 2374','Fe II 2382','Fe II 2586','Fe II 2600']}
    geometry=[dict(z=float(z),age_Gyr=float(age_gyr(z)),H2=float(h(z)),
                   example_relative_D_m_s=float(example_B*h(z)),
                   line_positions_nm={k:v*(1+z)/10 for k,v in wavelengths.items()}) for z in Z]
    # Independent numerical checks: Friedmann integral, affine prediction,
    # correlated-error cancellation and analytic Gaussian power threshold.
    checks=[]
    def check(name,condition):checks.append(dict(name=name,passed=bool(condition)))
    for z in [0.,*Z,3.]:
        amax=1/(1+z)
        integrated=quad(lambda a: np.sqrt(a)/np.sqrt(OM+(1-OM)*a**3),0,amax,
                        epsabs=1e-12,epsrel=1e-12)[0]*MPC_KM/H0/(YEAR*1e9)
        check(f'age_integral_z{z}',np.isclose(integrated,age_gyr(z),rtol=1e-10))
    check('closure_cancels_intercept',abs(contrast.sum())<1e-14)
    check('closure_cancels_p2_amplitude',abs(contrast@h(Z))<1e-14)
    check('constant_null_also_satisfies_closure',abs(contrast@np.ones(3))<1e-14)
    for shared_sigma in [0.,20.,100.]:
        cov=np.eye(3)*30**2+shared_sigma**2*np.ones((3,3))
        check(f'shared_lab_error_cancels_{shared_sigma}',
              np.isclose(contrast@cov@contrast,30**2*(contrast@contrast)))
    check('exact_gaussian_power_threshold',abs(power(threshold,1.)-.9)<1e-12)
    check('null_size',abs(power(0.,1.)-.01)<1e-14)
    check('no_advertised_n_below_power',all(
        s['required_n_per_endpoint_for_90pct_power_at_B100'] is None or
        power(example_B,np.sqrt(2*(s['single_target_random_sigma_m_s']**2/
              s['required_n_per_endpoint_for_90pct_power_at_B100']+
              s['nonaveraging_group_sigma_m_s']**2)))>=.9 for s in scenarios))
    result=dict(
        status='DESIGN ONLY; NOT PREREGISTERED; NO OBSERVED AMPLITUDE OR NEW DATA',
        model='D(z)=b+B H2(z); b=D(z=1), B=D(z=2)-D(z=1)',
        observable='D=v(Fe II2382)-v(Fe II2374); an apparent pair descriptor',
        physical_bridge='A deterministic universal mean response is an additional hypothesis; not derived from a Riemann resource/uncertainty floor.',
        cosmology=dict(H0_km_s_Mpc=H0,Omega_m=OM,flat=True,radiation=False,tstar_seconds=TSTAR),
        geometry=geometry,
        conditional_prediction=dict(
            anchor_z=[1.,2.],held_out_z=1.5,
            anchor_weights=[float(1-ratio),ratio],
            contrast_coefficients_low_mid_high=contrast.tolist(),
            variance='w.T @ full_joint_covariance @ w; include development, validation and their shared terms',
            independent_equal_group_sigma_multiplier=float(np.linalg.norm(contrast)),
            maximum_closure_bias_for_per_group_bound_epsilon=float(np.abs(contrast).sum()),
            example_endpoint_difference_m_s=example_B,
            example_middle_minus_low_m_s=example_B*ratio,
            example_is_not_measured_or_theory_predicted=True),
        competitors=alternatives,
        forecasts=dict(alpha_two_sided=.01,power_target=.9,
                       exact_required_effect_over_sigma=threshold,scenarios=scenarios,
                       caveat='Exact known-Gaussian planning model only. Group floors are hypothetical random offsets, not measured calibration priors. Unknown age-correlated biases need bounds or external measurements. Not a 5-sigma discovery calibration.'),
        operational_status=dict(target_ids_selected=False,exposure_time_calculated=False,
                                telescope_time_awarded=False,prediction_amplitude_frozen=False,
                                protocol_publicly_registered=False),
        source_hashes={'atomic_identification_csv':hashlib.sha256(atomic_path.read_bytes()).hexdigest(),
                       'design_code':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        validation=dict(passed=all(c['passed'] for c in checks),checks=checks))
    (ROOT/'results/design.json').write_text(json.dumps(result,indent=2)+'\n')
    assert result['validation']['passed']
    figure(result)
    print(json.dumps({'middle_fraction':ratio,'example_middle_m_s':example_B*ratio,
                      'competitors':alternatives,'checks':len(checks)},indent=2))


def figure(result):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.4))
    z=np.linspace(1.,2.,250)
    for p,style,color in [(1,'--','#c77922'),(3,':','#498355'),(2,'-','#245db1')]:
        axs[0].plot(z,100*h(z,p),style,color=color,label=f'p={p}',lw=1.8)
    mid=result['conditional_prediction']['example_middle_minus_low_m_s']
    axs[0].scatter([1.,2.],[0.,100.],marker='s',facecolors='white',edgecolors='black',s=55,zorder=5,
                   label='Illustrative anchors, not data')
    axs[0].scatter([1.5],[mid],color='#a73738',marker='D',s=40,zorder=5,label='Held-out prediction')
    axs[0].annotate(f'{mid:.2f} m/s',xy=(1.5,mid),xytext=(1.57,mid-18),
                    arrowprops=dict(arrowstyle='->',color='#555555'))
    axs[0].set(xlabel='Absorber redshift',ylabel='D(z) - D(1)  [m/s]',
               title='Conditional prediction: assumed B = 100 m/s')
    axs[0].legend(fontsize=8,loc='upper left')
    n=np.arange(2,61)
    for scenario,color in zip(result['forecasts']['scenarios'],['#245db1','#c77922','#498355']):
        s=scenario['single_target_random_sigma_m_s'];tau=scenario['nonaveraging_group_sigma_m_s']
        se=np.sqrt(2*(s*s/n+tau*tau))
        axs[1].plot(n,power(100,se),color=color,
                    label=rf'$\sigma_{{one}}={s:g},\ \tau_{{group}}={tau:g}$ m/s')
    axs[1].axhline(.9,color='gray',ls='--',lw=1,label='90% power')
    axs[1].set(xlabel='Independent targets per endpoint group',ylabel='Conditional detection power',
               ylim=(0,1.03),title='Planning only: B = 100 m/s, two-sided alpha = 0.01')
    axs[1].legend(fontsize=8,loc='lower right')
    fig.text(.5,.005,'No observations are plotted. Instrument calibration and gas-model adequacy are additional requirements.',
             ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.035,1,1])
    for ext in ['png','pdf']:fig.savefig(ROOT/f'figures/prediction_design.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':generate()
