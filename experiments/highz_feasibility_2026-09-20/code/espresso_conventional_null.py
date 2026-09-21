#!/usr/bin/env python3
"""Joint six-FeII conventional Voigt null and free relative-shift alternative.

A scoped independent analysis, NOT an exact reproduction of VPFIT or a cosmic
clock constraint. Fixed laboratory isotope wavelengths, abundance-weighted f,
natural damping, shared gas components, Gaussian instrumental convolution and
pixel integration. Data/model selection comes from the published right-region
model. Both hypotheses re-fit the same gas and continuum nuisance parameters.
"""
from pathlib import Path
import sys, re, json, argparse, hashlib, time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'.deps'))
import numpy as np
from astropy.io import fits
from scipy.special import wofz
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares
from scipy.stats import chi2
RAW=ROOT/'data/raw/espresso_null'; OUT=ROOT/'results/espresso_null'
C=299792.458; ZREF=1.150793
TAU_CONSTANT=1.497364150e-15  # sqrt(pi)*e^2/(m_e*c), with Angstrom, km/s, cm^-2
LINES=[2260,2344,2374,2382,2586,2600]; ANCHOR=2374
AUTHOR_ZERO={2260:0.,2344:.002,2374:0.,2382:.002,2586:0.,2600:.003}
AUTHOR_CONT={2260:1.,2344:1.,2374:.994,2382:1.,2586:.996,2600:1.}

def number(s):return float(re.match(r'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?',s).group())
def load_sources():
    s=(RAW/'fit_r_iso.f18').read_text(); last=s[s.rfind(' iteration   :'):].split('Parameter errors:')[0]
    comp=np.array([[number(x) for x in line.split()[1:4]] for line in last.splitlines() if line.strip().startswith('FeII ')])
    assert comp.shape==(45,3),comp.shape
    comp[:,1]=C*np.log((1+comp[:,1])/(1+ZREF))
    atoms={k:[] for k in LINES}
    for line in (RAW/'MM_VPFIT_2013-11-10.dat').read_text().splitlines():
        x=line.split()
        if len(x)>5 and x[0]=='FeII':
            for k in LINES:
                if k<=float(x[1])<k+1:atoms[k].append([float(t) for t in x[1:6]])
    atoms={k:np.array(v) for k,v in atoms.items()}
    regions={}
    for line in (RAW/'fit_r_iso.f13').read_text().splitlines():
        if line.startswith('hes'):
            for k in LINES:
                if f'FeII{k}' in line:
                    x=line.split();regions[k]=(float(x[2]),float(x[3]),float(x[4].split('=')[1]))
    with fits.open(RAW/'hes0515m4414.fits',memmap=False) as h:
        a=np.asarray(h[0].data,dtype=float);hdr=h[0].header
        wave=10**(hdr['CRVAL1']+(np.arange(a.shape[1])+1-hdr['CRPIX1'])*hdr['CD1_1'])
        dv=C*np.log(10)*hdr['CD1_1']
    return comp,atoms,regions,wave,a,dv

class Model:
    def __init__(self, free_shifts=False, continuum=True, error_row=1, oversample=9, fwhm_scale=1., rho=0., exclude=(), fit_velocity=True, fit_b=True, strength_nuisance=False, zero_nuisance=False, velocity_span=2., continuum_span=.03):
        comp,atoms,regions,wave,a,dv=load_sources()
        self.comp=comp;self.ncomp=len(comp);self.keys=[k for k in LINES if k not in exclude]
        self.free_shifts=free_shifts;self.continuum=continuum;self.fit_velocity=fit_velocity;self.fit_b=fit_b;self.rho=rho
        self.velocity_span=velocity_span;self.continuum_span=continuum_span
        self.labels=[f'logN_{i}' for i in range(len(comp))]
        if fit_velocity:self.labels += [f'v_{i}' for i in range(len(comp))]
        if fit_b:self.labels += [f'logb_{i}' for i in range(len(comp))]
        self.ngas=len(self.labels);self.ncont=2*len(self.keys) if continuum else 0
        self.strength_keys=[k for k in self.keys if k!=ANCHOR] if strength_nuisance else []
        self.zero_keys=self.keys.copy() if zero_nuisance else []
        self.strength_offset=self.ngas+self.ncont
        self.zero_offset=self.strength_offset+len(self.strength_keys)
        self.shift_offset=self.zero_offset+len(self.zero_keys)
        self.shift_keys=[k for k in self.keys if k!=ANCHOR] if free_shifts else []
        if free_shifts and ANCHOR not in self.keys:raise ValueError('The reference transition must remain in the analysis')
        self.labels += [f'{name}_{k}' for k in self.keys for name in ['continuum0','continuum1']] if continuum else []
        self.labels += [f'log_strength_{k}' for k in self.strength_keys]
        self.labels += [f'zero_adjustment_{k}' for k in self.zero_keys]
        self.labels += [f'shift_{k}' for k in self.shift_keys]
        self.npar=len(self.labels);self.lines=[];self.oversample=oversample;self.dv=dv
        self.slices=[];n=0
        for ki,k in enumerate(self.keys):
            atom=atoms[k];ref=np.sum(atom[:,0]*atom[:,1])/atom[:,1].sum()
            lo,hi,fwhm=regions[k]
            use=(wave>=lo)&(wave<=hi)
            idx=np.flatnonzero(use)
            good=(a[4,idx]==1)&(a[error_row,idx]>0)&np.isfinite(a[0,idx])
            v=C*np.log(wave[idx]/(ref*(1+ZREF)))
            pad=int(np.ceil(6*fwhm/dv))+2
            grid=(v[0]-pad*dv)+(np.arange((len(idx)+2*pad)*oversample)+.5)*dv/oversample-.5*dv
            # Subpixel means exactly centered on data pixel centers.
            segments=np.flatnonzero(np.r_[True,np.diff(idx[good])>1])
            self.lines.append(dict(key=k,atom=atom,ref=ref,v=v,grid=grid,pad=pad,good=good,
                flux=a[0,idx],error=a[error_row,idx],idx=idx,fwhm=fwhm*fwhm_scale,segments=segments,
                x=(v-v.mean())/100.,wave=wave[idx]))
            self.slices.append(slice(n,n+good.sum()));n+=good.sum()
        self.ndata=int(n);self.cache=None
    def initial(self):
        c=self.comp
        parts=[c[:,0]];lo=[np.maximum(7.,c[:,0]-2.)];hi=[np.minimum(17.,c[:,0]+2.)]
        if self.fit_velocity:
            parts+=[c[:,1]];lo+=[c[:,1]-self.velocity_span];hi+=[c[:,1]+self.velocity_span]
        if self.fit_b:
            parts+=[np.log(c[:,2])];lo+=[np.log(np.full(len(c),.15))];hi+=[np.log(np.full(len(c),30.))]
        if self.continuum:parts+=[np.zeros(self.ncont)];lo+=[np.full(self.ncont,-self.continuum_span)];hi+=[np.full(self.ncont,self.continuum_span)]
        if self.strength_keys:parts+=[np.zeros(len(self.strength_keys))];lo+=[np.full(len(self.strength_keys),np.log(.9))];hi+=[np.full(len(self.strength_keys),np.log(1.1))]
        if self.zero_keys:parts+=[np.zeros(len(self.zero_keys))];lo+=[np.full(len(self.zero_keys),-.01)];hi+=[np.full(len(self.zero_keys),.01)]
        if self.free_shifts:parts+=[np.zeros(len(self.shift_keys))];lo+=[np.full(len(self.shift_keys),-.5)];hi+=[np.full(len(self.shift_keys),.5)]
        return np.concatenate(parts),np.concatenate(lo),np.concatenate(hi)
    def unpack(self,p):
        n=self.ncomp;j=n;N=10**p[:n]
        v=p[j:j+n] if self.fit_velocity else self.comp[:,1]
        j+=n if self.fit_velocity else 0
        b=np.exp(p[j:j+n]) if self.fit_b else self.comp[:,2]
        return N,v,b
    def evaluate(self,p,jac=True):
        if self.cache is not None and np.array_equal(p,self.cache[0]):return self.cache[1:]
        N,vc,b=self.unpack(p);n=self.ncomp;res=[];jacrows=[];models=[]
        for ki,L in enumerate(self.lines):
            k=L['key'];g=L['grid'];shift=p[self.shift_offset+self.shift_keys.index(k)] if k in self.shift_keys else 0.
            strength=np.exp(p[self.strength_offset+self.strength_keys.index(k)]) if k in self.strength_keys else 1.
            tau=np.zeros((len(g),n));tv=np.zeros_like(tau);tb=np.zeros_like(tau)
            for lam,f,gamma,mass,q in L['atom']:
                ratio=np.exp((vc[None,:]+C*np.log(lam/L['ref'])+shift-g[:,None])/C)
                u=C*(1-ratio)/b[None,:]
                damping=gamma*lam*1e-13/(4*np.pi*b)
                z=u+1j*damping[None,:];w=wofz(z);wp=-2*z*w+2j/np.sqrt(np.pi)
                amp=TAU_CONSTANT*N*f*lam/b*strength
                t=w.real*amp[None,:];tau+=t
                tv+=amp[None,:]*wp.real*(-ratio/b[None,:])
                tb+=-t+amp[None,:]*(-u*wp.real+damping[None,:]*wp.imag)
            intrinsic=np.exp(-tau.sum(axis=1))
            cols=[-intrinsic[:,None]*tau*np.log(10)]
            if self.fit_velocity:cols+=[-intrinsic[:,None]*tv]
            if self.fit_b:cols+=[-intrinsic[:,None]*tb]
            cols+=[(-intrinsic*tv.sum(axis=1))[:,None]]
            cols+=[(-intrinsic*tau.sum(axis=1))[:,None]]
            allarr=np.column_stack([intrinsic,*cols])
            conv=gaussian_filter1d(allarr,L['fwhm']/2.354820045/(self.dv/self.oversample),axis=0,mode='nearest',truncate=6.)
            pix=conv.reshape(-1,self.oversample,conv.shape[1]).mean(axis=1)[L['pad']:-L['pad']]
            c0,c1=p[self.ngas+2*ki:self.ngas+2*ki+2] if self.continuum else (0.,0.)
            cont=AUTHOR_CONT[k]*(1+c0+c1*L['x'])
            zero=AUTHOR_ZERO[k]+(p[self.zero_offset+self.zero_keys.index(k)] if k in self.zero_keys else 0.)
            profile=zero+(1-zero)*pix[:,0]
            m=cont*profile;models.append(m)
            J=np.zeros((len(m),self.npar));J[:,:self.ngas]=(cont*(1-zero))[:,None]*pix[:,1:1+self.ngas]
            if self.continuum:
                J[:,self.ngas+2*ki]=AUTHOR_CONT[k]*profile;J[:,self.ngas+2*ki+1]=AUTHOR_CONT[k]*L['x']*profile
            if k in self.strength_keys:J[:,self.strength_offset+self.strength_keys.index(k)]=cont*(1-zero)*pix[:,-1]
            if k in self.zero_keys:J[:,self.zero_offset+self.zero_keys.index(k)]=cont*(1-pix[:,0])
            if k in self.shift_keys:J[:,self.shift_offset+self.shift_keys.index(k)]=cont*(1-zero)*pix[:,-2]
            good=L['good'];r=(m[good]-L['flux'][good])/L['error'][good];J=J[good]/L['error'][good,None]
            if self.rho:
                # Fixed AR(1) sensitivity model, reset after every masked gap.
                previous=np.r_[False,np.diff(L['idx'][good])==1]
                r0=r.copy();J0=J.copy();r[previous]=(r0[previous]-self.rho*r0[np.flatnonzero(previous)-1])/np.sqrt(1-self.rho**2)
                J[previous]=(J0[previous]-self.rho*J0[np.flatnonzero(previous)-1])/np.sqrt(1-self.rho**2)
            res.append(r);jacrows.append(J)
        self.cache=(p.copy(),np.concatenate(res),np.vstack(jacrows),models)
        return self.cache[1:]
    def fun(self,p):
        r=self.evaluate(p)[0]
        if hasattr(self,'progress_name'):
            self.progress_count+=1
            if self.progress_count%25==0:
                print(f'{self.progress_name}: evaluation {self.progress_count}, chi2={r@r:.7f}',flush=True)
                checkpoint=dict(labels=self.labels,parameters=p.tolist(),chi2=float(r@r),evaluation=self.progress_count)
                (OUT/f'{self.progress_name}_checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
        return r
    def jac(self,p):return self.evaluate(p)[1]

def fit_case(name,free=False,start=None,max_nfev=300,**kwargs):
    OUT.mkdir(parents=True,exist_ok=True);model=Model(free_shifts=free,**kwargs)
    p,lo,hi=model.initial()
    if start is not None:
        old=json.loads(Path(start).read_text());lookup=dict(zip(old['labels'],old['parameters']))
        p=np.array([lookup.get(k,val) for k,val in zip(model.labels,p)])
    p=np.clip(p,lo+1e-9,hi-1e-9);t=time.time();r0=model.fun(p)
    print(f'{name}: ndata={model.ndata} npar={model.npar} initial chi2={r0@r0:.6f}',flush=True)
    model.progress_name=name;model.progress_count=0
    opt=least_squares(model.fun,p,jac=model.jac,bounds=(lo,hi),x_scale='jac',ftol=1e-7,xtol=1e-8,gtol=1e-6,max_nfev=max_nfev,verbose=0)
    r,J,profiles=model.evaluate(opt.x);sv=np.linalg.svd(J,compute_uv=False);rank=int(sum(sv>sv[0]*1e-10))
    _,ss,Vh=np.linalg.svd(J,full_matrices=False)
    keep=ss>ss[0]*1e-10
    covariance=(Vh[keep].T/ss[keep]**2)@Vh[keep]
    # The direct-SVD covariance is a local diagnostic only; poorly identified directions and
    # bound-active nuisance parameters make naive Gaussian inference unsafe.
    chi=float(r@r);ndf=model.ndata-model.npar
    out=dict(name=name,free_shifts=free,configuration=kwargs,ndata=model.ndata,npar=model.npar,nominal_ndf=ndf,chi2=chi,
        chi2_per_nominal_ndf=chi/ndf,nominal_chi2_survival=float(chi2.sf(chi,ndf)),
        optimizer_success=bool(opt.success),message=opt.message,nfev=opt.nfev,optimality=float(opt.optimality),seconds=time.time()-t,
        jacobian_rank=rank,jacobian_singular_values=sv.tolist(),labels=model.labels,parameters=opt.x.tolist(),
        active_bounds=[model.labels[i] for i in range(model.npar) if min(opt.x[i]-lo[i],hi[i]-opt.x[i])<1e-5],
        shifts_km_s={k:float(opt.x[model.labels.index(f'shift_{k}')]) for k in model.shift_keys},
        per_line=[dict(line=L['key'],ndata=int(L['good'].sum()),chi2=float(r[sl]@r[sl]),chi2_per_pixel=float(r[sl]@r[sl]/L['good'].sum())) for L,sl in zip(model.lines,model.slices)],
        interpretation='Conditional six-FeII, one-absorber test; no cosmic-time law fitted. Nominal diagonal-noise metrics are not discovery significances.')
    (OUT/f'{name}.json').write_text(json.dumps(out,indent=2)+'\n')
    arrays=dict(parameters=opt.x,residuals=r,jacobian=J,covariance_local_diagnostic=covariance)
    for L,m in zip(model.lines,profiles):
        for key in ['v','wave','flux','error','good']:arrays[f'{L["key"]}_{key}']=L[key]
        arrays[f'{L["key"]}_model']=m
    np.savez_compressed(OUT/f'{name}.npz',**arrays)
    print(json.dumps({k:v for k,v in out.items() if k not in ['labels','parameters','jacobian_singular_values']},indent=2),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--name',default='null');p.add_argument('--free',action='store_true');p.add_argument('--start');p.add_argument('--rho',type=float,default=0.);p.add_argument('--fwhm-scale',type=float,default=1.);p.add_argument('--error-row',type=int,default=1);p.add_argument('--max-nfev',type=int,default=300);p.add_argument('--oversample',type=int,default=9);p.add_argument('--exclude',type=int,nargs='*',default=[]);p.add_argument('--fixed-continuum',action='store_true');p.add_argument('--strength-nuisance',action='store_true');p.add_argument('--zero-nuisance',action='store_true')
    a=p.parse_args();fit_case(a.name,a.free,a.start,max_nfev=a.max_nfev,rho=a.rho,fwhm_scale=a.fwhm_scale,error_row=a.error_row,oversample=a.oversample,exclude=a.exclude,continuum=not a.fixed_continuum,strength_nuisance=a.strength_nuisance,zero_nuisance=a.zero_nuisance)
