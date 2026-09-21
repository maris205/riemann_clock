#!/usr/bin/env python3
"""Conditional UVES SQUAD cross-instrument fit using isotope Voigt physics.

This has not been established as the final Kotus+2017 UVES product. The UVES coadd has no measured LSF
or complete wavelength-distortion/covariance product. Relative shifts are a
replication diagnostic conditional on this forward model, not a cosmic result.
The baseline excludes known-blended 2344/2586 and weak 2260; optional inclusion
of 2586 uses author ESPRESSO mask wavelengths expanded by 3 km/s as a conservative
shared-blend control. The extra expansion is a declared sensitivity choice,
not an assertion that ESPRESSO instrumental bad pixels transfer to UVES.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import json, argparse, hashlib
import numpy as np
import espresso_conventional_null as base
ROOT=base.ROOT
OUT=ROOT/'results/cross_instrument'
ORIGINAL_LOAD=base.load_sources
ORIGINAL_MODEL=base.Model
ORIGINAL_FILTER=base.gaussian_filter1d

def load_uves(fwhm_scale=1., transfer_masks=True):
    comp,atoms,regions,ewave,ea,edv=ORIGINAL_LOAD()
    z=np.load(ROOT/'data/processed/J051707-441055_squad_dr1.npz')
    w=z['wavelength_vacuum_heliocentric_AA']
    a=np.zeros((5,len(w)))
    a[0]=z['flux_normalized'];a[1]=z['error_normalized'];a[2]=z['expected_fluctuation_normalized'];a[3]=1.;a[4]=z['valid']
    dv=float(np.median(base.C*np.diff(np.log(w))))
    masks=json.loads((base.RAW/'feii_mask_audit.json').read_text())
    removed=[]
    if transfer_masks:
        for row in masks:
            if '2586' not in row['region']['description']:continue
            for span in row['masked_spans']:
                lo=span['w_start_AA']*np.exp(-(3.+edv/2)/base.C)
                hi=span['w_end_AA']*np.exp((3.+edv/2)/base.C)
                bad=(w>=lo)&(w<=hi);a[4,bad]=0
                removed.append(dict(line=2586,start_AA=lo,end_AA=hi,uv_pixels=int(bad.sum())))
    # SQUAD nominal resolution in this wavelength band, not a measured LSF.
    regions={k:(v[0],v[1],base.C/53696.0*fwhm_scale) for k,v in regions.items()}
    OUT.mkdir(parents=True,exist_ok=True)
    info=dict(instrument='UVES',source_product='UVES SQUAD DR1 public coadd; exact equivalence to final Kotus+2017 corrected spectrum has not been verified',
              nominal_R=53696.,fwhm_km_s=base.C/53696.*fwhm_scale,pixel_km_s=dv,
              relative_mask_origin='ESPRESSO published masked FeII2586 spans expanded by 3 km/s plus half ESPRESSO pixel',transferred_mask_spans=removed,
              source_npz_sha256=hashlib.sha256((ROOT/'data/processed/J051707-441055_squad_dr1.npz').read_bytes()).hexdigest(),
              historical_observation_dates=['1999-12-14','2009-02-11'],wavelength='vacuum heliocentric Angstrom')
    (OUT/'uves_data_metadata.json').write_text(json.dumps(info,indent=2)+'\n')
    return comp,atoms,regions,w,a,dv

class FlexibleLSFModel(ORIGINAL_MODEL):
    """Same gas, one Gaussian LSF width per transition as conventional nuisance."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.base_npar=self.npar
        self.labels += [f'log_lsf_scale_{k}' for k in self.keys]
        self.npar=len(self.labels)
        self.nominal_fwhm=[L['fwhm'] for L in self.lines]
        self.extended_cache=None
    def initial(self):
        p,lo,hi=super().initial()
        return np.r_[p,np.zeros(len(self.keys))],np.r_[lo,np.full(len(self.keys),np.log(.65))],np.r_[hi,np.full(len(self.keys),np.log(1.25))]
    def evaluate(self,p,jac=True):
        if self.extended_cache is not None and np.array_equal(p,self.extended_cache[0]):return self.extended_cache[1:]
        for i,L in enumerate(self.lines):L['fwhm']=self.nominal_fwhm[i]*np.exp(p[self.base_npar+i])
        # Capture only the intrinsic-flux convolution derivative. This gives the
        # same centered finite difference as re-running the full Voigt calculation,
        # while avoiding six duplicate optical-depth/Jacobian calculations.
        width_derivatives=[]
        step=2e-5
        def capture(arr,sigma,**kw):
            out=ORIGINAL_FILTER(arr,sigma,**kw)
            plus=ORIGINAL_FILTER(arr[:,0],sigma*np.exp(step),**kw)
            minus=ORIGINAL_FILTER(arr[:,0],sigma*np.exp(-step),**kw)
            width_derivatives.append((plus-minus)/(2*step))
            return out
        self.cache=None
        base.gaussian_filter1d=capture
        try:
            r,J,profiles=super().evaluate(p)
        finally:
            base.gaussian_filter1d=ORIGINAL_FILTER
        r=r.copy();J=J.copy();profiles=[m.copy() for m in profiles]
        for i,(L,dconv,sl) in enumerate(zip(self.lines,width_derivatives,self.slices)):
            dpix=dconv.reshape(-1,self.oversample).mean(axis=1)[L['pad']:-L['pad']]
            c0,c1=p[self.ngas+2*i:self.ngas+2*i+2] if self.continuum else (0.,0.)
            cont=base.AUTHOR_CONT[L['key']]*(1+c0+c1*L['x'])
            zero=base.AUTHOR_ZERO[L['key']]+(p[self.zero_offset+self.zero_keys.index(L['key'])] if L['key'] in self.zero_keys else 0.)
            col=(cont*(1-zero)*dpix)[L['good']]/L['error'][L['good']]
            if self.rho:
                previous=np.r_[False,np.diff(L['idx'][L['good']])==1]
                original=col.copy()
                col[previous]=(original[previous]-self.rho*original[np.flatnonzero(previous)-1])/np.sqrt(1-self.rho**2)
            J[sl,self.base_npar+i]=col

        self.cache=None
        self.extended_cache=(p.copy(),r,J,profiles)
        return self.extended_cache[1:]

def configure(fwhm_scale=1.):
    data=load_uves(fwhm_scale)
    base.load_sources=lambda:data
    base.OUT=OUT
    base.AUTHOR_CONT={k:1. for k in base.LINES}
    base.AUTHOR_ZERO={k:0. for k in base.LINES}

def run(name,free=False,start=None,include_blended=False,fwhm_scale=1.,rho=0.,max_nfev=180,flex_lsf=False,oversample=21,error_row=1,strength_nuisance=False):
    configure(fwhm_scale)
    base.Model=FlexibleLSFModel if flex_lsf else ORIGINAL_MODEL
    out=base.fit_case(name,free=free,start=start,max_nfev=max_nfev,oversample=oversample,
        exclude=[2260,2344] if include_blended else [2260,2344,2586],rho=rho,error_row=error_row,
        zero_nuisance=True,strength_nuisance=strength_nuisance,velocity_span=5.,continuum_span=.03)
    out['cross_instrument_configuration']=dict(flexible_lsf=flex_lsf,nominal_fwhm_scale=fwhm_scale,include_blended=include_blended)
    out['interpretation']='Conditional UVES SQUAD DR1 cross-instrument diagnostic; catalogue LSF, incompletely audited wavelength-correction provenance, and shared gas model limit replication. No cosmic-time law is fitted.'
    out['input_start_file']=str(start) if start else None
    out['source_code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (OUT/f'{name}.json').write_text(json.dumps(out,indent=2)+'\n')
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--name',default='uves_clean_null');p.add_argument('--free',action='store_true');p.add_argument('--start');p.add_argument('--include-blended',action='store_true');p.add_argument('--fwhm-scale',type=float,default=1.);p.add_argument('--rho',type=float,default=0.);p.add_argument('--max-nfev',type=int,default=180);p.add_argument('--flex-lsf',action='store_true');p.add_argument('--oversample',type=int,default=21);p.add_argument('--error-row',type=int,default=1);p.add_argument('--strength-nuisance',action='store_true')
    a=p.parse_args();run(a.name,a.free,a.start,a.include_blended,a.fwhm_scale,a.rho,a.max_nfev,a.flex_lsf,a.oversample,a.error_row,a.strength_nuisance)
