#!/usr/bin/env python3
"""Exploratory registration of ACTUAL Fe II profiles, not a physical drift fit.

Use Fe II 2374 as a noisy empirical template for 2586 and 2600 in the same
archived spectrum. Power/affine transforms are descriptive and generally do not
commute with the LSF convolution. Neither registration shifts nor formal errors
are alpha, cosmic-aging or Riemann-cutoff estimates.
"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.optimize import minimize, differential_evolution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
TARGET = 'J051707-441055'
ZABS = 1.1508
C_KMS = 299792.458
FILE_LABELS = {2374:'2374.4601',2586:'2586.6493',2600:'2600.1725'}
REGIONS = {'whole':(-25.,110.), 'central':(-15.,40.), 'red':(40.,95.)}
ERRORS = {'statistical':'error_normalized','fluctuation':'expected_fluctuation_normalized'}
BOUNDS = [(-3.,3.),(.15,8.),(.90,1.10),(-.15,.15)]


def read_data():
    lines={int(r['label'].split()[-1]):float(r['ritz_vacuum_wavelength_A'])
           for r in csv.DictReader((BASE/'data/atomic/selected_transitions.csv').open()) if r['ion']=='Fe II'}
    arrays={}; hashes={}
    for label,filename in FILE_LABELS.items():
        path=BASE/'data/processed'/f'{TARGET}_FeII_{filename}_window.npz'
        hashes[str(path.relative_to(BASE))]=hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path,allow_pickle=False) as npz:
            data={k:npz[k] for k in npz.files}
        data['v']=C_KMS*(data['wavelength_vacuum_heliocentric_AA']/(lines[label]*(1+ZABS))-1)
        arrays[label]=data
    return lines, arrays, hashes


def interp_template(v,ref,error):
    i=np.searchsorted(ref['v'],v)-1
    i=np.clip(i,0,len(ref['v'])-2)
    w=(v-ref['v'][i])/(ref['v'][i+1]-ref['v'][i])
    if np.any((w<0)|(w>1)):
        raise ValueError('Template extrapolation')
    f=(1-w)*ref['flux_normalized'][i]+w*ref['flux_normalized'][i+1]
    # Diagonal approximation only. Interpolation also induces cross-row
    # covariance which is NOT supplied by this registration likelihood.
    var=(1-w)**2*ref[error][i]**2+w**2*ref[error][i+1]**2
    return f,var


def prediction(p,v,ref,error,form):
    shift, strength, c0, c1=p
    f,var=interp_template(v-shift,ref,error)
    if np.any(f<=0):
        raise ValueError('Nonpositive template flux in selected domain')
    continuum=c0+c1*v/100.
    if form=='power':
        model=continuum*f**strength
        derivative=continuum*strength*f**(strength-1)
    else:
        model=continuum+strength*(f-1)
        derivative=np.full_like(f,strength)
    return model,derivative**2*var


def finite_hessian(fun,p):
    steps=np.array([.002,.0003,.00003,.00003])
    h=np.zeros((4,4)); f0=fun(p)
    for i in range(4):
        di=np.eye(4)[i]*steps[i]
        h[i,i]=(fun(p+di)-2*f0+fun(p-di))/steps[i]**2
        for j in range(i):
            dj=np.eye(4)[j]*steps[j]
            h[i,j]=h[j,i]=(fun(p+di+dj)-fun(p+di-dj)-fun(p-di+dj)+fun(p-di-dj))/(4*steps[i]*steps[j])
    return h


def fit_case(ref,data,region,error_name,form,label):
    low,high=REGIONS[region]; error=ERRORS[error_name]
    valid=data['valid']&(data['v']>=low)&(data['v']<=high)&(data[error]>0)
    v=data['v'][valid]; y=data['flux_normalized'][valid]; s=data[error][valid]
    def objective(p):
        m,rv=prediction(p,v,ref,error,form)
        variance=s*s+rv
        return float(np.sum((y-m)**2/variance+np.log(variance)))
    fits=[]
    for start_shift in (-.5,0.,.5):
        p0=np.array([start_shift,2. if label==2586 else 5.,1.,0.])
        fit=minimize(objective,p0,method='Powell',bounds=BOUNDS,
                     options={'xtol':1e-8,'ftol':1e-10,'maxiter':600})
        fits.append(fit)
    best=min(fits,key=lambda f:f.fun)
    start_range=float(max(f.fun for f in fits)-min(f.fun for f in fits))
    # A local fit can settle at different interpolation cells. Check every
    # configuration globally, including cases whose local starts agree.
    global_check=True
    global_values=[]
    for seed in (726308,25862600):
        de=differential_evolution(objective,BOUNDS,seed=seed,popsize=20,maxiter=500,tol=1e-10,polish=True)
        refined=minimize(objective,de.x,method='Powell',bounds=BOUNDS,
                         options={'xtol':1e-9,'ftol':1e-11,'maxiter':800})
        global_values.append(float(min(de.fun,refined.fun)))
        best=min([best,de,refined],key=lambda f:f.fun)
    p=best.x; m,rv=prediction(p,v,ref,error,form); var=s*s+rv
    h=finite_hessian(objective,p)
    eig=np.linalg.eigvalsh(h)
    boundary=any(min(x-lo,hi-x)<1e-4*(hi-lo) for x,(lo,hi) in zip(p,BOUNDS))
    # Linear interpolation has derivative discontinuities. A Hessian evaluated
    # across a knot can manufacture a tiny "error"; withhold it in that case.
    shifted=v-p[0]
    nearest=np.min(np.abs(shifted[:,None]-ref['v'][None,:]),axis=1)
    knot_distance=float(nearest.min())
    smooth=knot_distance>0.004
    sigma=float(np.sqrt(2*np.linalg.inv(h)[0,0])) if np.all(eig>0) and smooth and not boundary and best.success else None
    chi2=float(np.sum((y-m)**2/var)); dof=len(y)-len(p)
    summary={'target':TARGET,'reference_FeII':2374,'target_FeII':label,'region':region,
        'velocity_bounds_kms':[low,high],'error_array':error_name,'form':form,'n_pixels':len(v),
        'shift_kms':float(p[0]),'formal_sigma_shift_kms':sigma,
        'strength':float(p[1]),'continuum_intercept':float(p[2]),'continuum_slope_per100kms':float(p[3]),
        'chi2_diagonal':chi2,'nominal_dof':dof,'chi2_per_nominal_dof':chi2/dof,
        'gaussian_deviance_including_log_variance':float(best.fun),
        'optimizer_success':bool(best.success),'parameter_at_boundary':boundary,
        'positive_local_hessian':bool(np.all(eig>0)),
        'min_distance_to_template_interpolation_knot_kms':knot_distance,
        'smooth_for_local_hessian':smooth,
        'formal_error_status':'conditional_curvature_only' if sigma is not None else 'WITHHELD_BOUNDARY_KNOT_OR_HESSIAN',
        'start_solution_deviance_range':start_range,
        'additional_global_search':global_check,
        'global_seed_deviance_range':float(np.ptp(global_values)),
        'interpretation':'DESCRIPTIVE_REGISTRATION_ONLY; not a physical frequency shift; independent-pixel approximation'}
    return summary,(v,y,s,m,np.sqrt(var))


def run():
    lines,arrays,hashes=read_data()
    rows=[]; profiles={}
    for label in (2586,2600):
        for region in REGIONS:
            for error_name in ERRORS:
                for form in ('power','affine'):
                    row,curve=fit_case(arrays[2374],arrays[label],region,error_name,form,label)
                    rows.append(row)
                    if region=='whole' and error_name=='fluctuation' and form=='power':
                        profiles[label]=curve
    summary={'status':'ACTUAL_OBSERVED_FLUX; exploratory uncorrected profile registration, no cosmic fit',
        'target':TARGET,'z_abs_reference':ZABS,'reference_transition':2374,
        'rest_vacuum_wavelength_A':lines,'source_window_sha256':hashes,
        'known_limitations':['Power of convolved flux is not a convolved optical-depth model',
            'Fe II 2586 has published additional absorption; this is an uncorrected stress test',
            'Fe II 2600 is strongly saturated; affine and power forms need not fit it',
            'Pixel and interpolation correlations are omitted; errors are local formal diagnostics only',
            'NIST laboratory and wavelength-calibration covariance omitted',
            'Regions were selected after inspecting the data, not preregistered or independent trials',
            'One local section of an extended absorber, not its complete ~720 km/s profile'],
        'fits':rows}
    (BASE/'results/observed_registration.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    csvrows=[{k:v for k,v in row.items() if k!='velocity_bounds_kms'} for row in rows]
    with (BASE/'results/observed_registration.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(csvrows[0]));w.writeheader();w.writerows(csvrows)
    fig,axes=plt.subplots(2,2,figsize=(11,7),gridspec_kw={'height_ratios':[2,1]},constrained_layout=True)
    for col,label in enumerate((2586,2600)):
        v,y,s,m,total=profiles[label]
        ax=axes[0,col]
        ax.errorbar(v,y,yerr=s,fmt='.',ms=3,lw=.6,color='#275f80',label='Observed target flux')
        ax.plot(v,m,color='#bc693e',lw=1.4,label='2374 template, power transform')
        row=next(r for r in rows if r['target_FeII']==label and r['region']=='whole' and r['error_array']=='fluctuation' and r['form']=='power')
        ax.set(title=f'Fe II {label}: descriptive registration',ylabel='Normalized observed flux')
        ax.text(.03,.07,f"Diagonal residual statistic / dof = {row['chi2_per_nominal_dof']:.1f}",transform=ax.transAxes,fontsize=9)
        ax.legend(fontsize=8,frameon=False)
        axes[1,col].axhline(0,color='grey',lw=.8)
        axes[1,col].plot(v,(y-m)/total,'.',color='#275f80',ms=3)
        axes[1,col].set(xlabel='Velocity relative to z=1.1508 (km/s)',ylabel='Residual / nominal error')
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle('Real HE 0515−4414 data: a rigid empirical template does not fit all line structure',fontsize=12)
    for suffix in ('png','pdf'):fig.savefig(BASE/f'figures/observed_registration.{suffix}',dpi=190)
    plt.close(fig)
    print(json.dumps([r for r in rows if r['error_array']=='fluctuation' and r['form']=='power'],indent=2))


if __name__=='__main__':
    run()
