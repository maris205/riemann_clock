#!/usr/bin/env python3
"""Summarize completed UVES fits without promoting nominal metrics to evidence."""
from pathlib import Path
import json, csv, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'results/cross_instrument'
def summarize():
    fits=[]
    for p in sorted(D.glob('uves_*.json')):
        if 'checkpoint' in p.name or 'metadata' in p.name:continue
        a=json.loads(p.read_text())
        if 'parameters' not in a or 'optimizer_success' not in a:continue
        vals=dict(zip(a['labels'],a['parameters']))
        row=dict(name=p.stem,free_shifts=a['free_shifts'],chi2=a['chi2'],ndata=a['ndata'],npar=a['npar'],rank=a['jacobian_rank'],
          success=a['optimizer_success'],nfev=a['nfev'],optimality=a['optimality'],oversample=a['configuration'].get('oversample'),
          error_row=a['configuration'].get('error_row',1),shift2382_m_s=1000*vals.get('shift_2382',0),
          shift2600_m_s=1000*vals.get('shift_2600',0),bound_count=len(a['active_bounds']),
          lsf_fwhm_km_s={k:299792.458/53696*np.exp(vals.get(f'log_lsf_scale_{k}',0)) for k in [2374,2382,2600]},
          source_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
        fits.append(row)
    out=dict(fits=fits,scope='Conditional independent-instrument data test with shared ESPRESSO-informed component architecture; not independent gas-model selection.',
             inference='No calibrated global significance, no alpha estimate, no inverse-log-squared time fit. Residual wavelength-calibration uncertainty and covariance not established.')
    for kind in ['refined','expected','opacity']:
        group=[f for f in fits if f['name'].startswith(f'uves_{kind}_') and f['oversample']>=21]
        nn=[f for f in group if not f['free_shifts']];aa=[f for f in group if f['free_shifts']]
        if nn and aa:
            n=min(nn,key=lambda f:f['chi2']);a=min(aa,key=lambda f:f['chi2'])
            out[kind]=dict(null=n['name'],alternative=a['name'],delta_chi2=n['chi2']-a['chi2'],extra_parameters=a['npar']-n['npar'],
                          both_optimizers_converged=n['success'] and a['success'],candidates=[f['name'] for f in group])
    (D/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
    if fits:
        with (D/'summary.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(fits[0]));w.writeheader();w.writerows(fits)
    p=D/'uves_refined_alternative.npz'
    if p.exists():
        z=np.load(p);fig,ax=plt.subplots(3,2,figsize=(11,8),gridspec_kw={'width_ratios':[2.1,1]})
        for i,k in enumerate([2374,2382,2600]):
            v=z[f'{k}_v'];fl=z[f'{k}_flux'];m=z[f'{k}_model'];er=z[f'{k}_error'];g=z[f'{k}_good']
            ax[i,0].step(v,fl,where='mid',color='#263d51',lw=.8,label='UVES SQUAD DR1')
            ax[i,0].plot(v,m,color='#dc6b2f',lw=1.4,label='H1: shared gas, LSF and relative shifts')
            ax[i,0].fill_between(v,m-er,m+er,color='#dc6b2f',alpha=.16)
            ax[i,0].set_ylabel(f'Fe II {k}\nNormalized flux');ax[i,0].set_ylim(-.08,1.13)
            ax[i,1].step(v[g],((fl-m)/er)[g],where='mid',color='#263d51',lw=.7)
            ax[i,1].axhline(0,color='gray',lw=.8);ax[i,1].axhspan(-1,1,color='gray',alpha=.15);ax[i,1].set_ylim(-5,5);ax[i,1].set_ylabel('Residual / supplied error')
        ax[0,0].legend(fontsize=8,loc='lower right');ax[-1,0].set_xlabel('Velocity relative to z = 1.150793 (km/s)');ax[-1,1].set_xlabel('Velocity (km/s)')
        fig.suptitle('Independent UVES data: a shared gas model with free relative line shifts',fontsize=12)
        fig.text(.5,.015,'Historical coadd; 45 gas components; H1 diagnostic fit; Gaussian widths are nuisance parameters.',ha='center',fontsize=9)
        fig.tight_layout(rect=(0,.04,1,.96));fig.savefig(D/'uves_profiles.png',dpi=180);fig.savefig(D/'uves_profiles.pdf');plt.close(fig)
    return out
if __name__=='__main__':print(json.dumps(summarize(),indent=2))
