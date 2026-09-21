#!/usr/bin/env python3
"""Exploratory physical-profile feasibility fits to two archival absorbers.
Shared Voigt opacity, terrestrial isotope mixture, variable per-line Gaussian
LSF, finite pixel integration, same fixed pixels for H0/H1. No alpha/time fit.
"""
from pathlib import Path
import json,hashlib,time,argparse
import numpy as np
from scipy.special import wofz
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares
from concurrent.futures import ProcessPoolExecutor
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/processed/archive_expansion'
OUT=ROOT/'results/archive_expansion/profile_pilot'
C=299792.458; TAU=1.497364150e-15
ATOMIC=ROOT/'data/raw/espresso_null/MM_VPFIT_2013-11-10.dat'
CONFIG={
'J232128-105122':dict(z=1.629,lines=[2260,2344,2374,2382,2586],window=[-30.,35.],component_bounds=[-20.,30.],centers=[2.,7.,12.,17.,22.,27.],logtotal=14.3),
'J004131-493611':dict(z=2.248,lines=[1608,2374,2382],window=[-30.,85.],component_bounds=[-25.,75.],centers=[-10.,20.,30.,40.,49.,57.],logtotal=14.4)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atoms(keys):
 d={k:[] for k in keys}
 for line in ATOMIC.read_text().splitlines():
  x=line.split()
  if len(x)>5 and x[0]=='FeII':
   for k in keys:
    if k<=float(x[1])<k+1:d[k].append([float(t) for t in x[1:6]])
 return {k:np.array(v) for k,v in d.items()}
class Model:
 def __init__(self,target,ncomp,free=False,oversample=11):
  self.target=target;self.cfg=CONFIG[target];self.ncomp=ncomp;self.free=free;self.oversample=oversample
  self.keys=self.cfg['lines'];self.shiftkeys=[k for k in self.keys if k!=2374] if free else []
  self.labels=[f'logN_{i}' for i in range(ncomp)]+[f'v_{i}' for i in range(ncomp)]+[f'logb_{i}' for i in range(ncomp)]
  self.labels += [f'{a}_{k}' for k in self.keys for a in ['cont0','cont1','zero','logfwhm']]
  self.labels += [f'shift_{k}' for k in self.shiftkeys]
  self.npar=len(self.labels);self.lines=[];self.slices=[];self.cache=None;off=0
  ad=atoms(self.keys);meta=json.loads((DATA/f'{target}_metadata.json').read_text());rows=meta['tables']['1']['rows']
  for k in self.keys:
   path=DATA/f'{target}_z{self.cfg["z"]:.3f}_FeII{k}_window.npz';d=np.load(path);a=ad[k];assert len(a)>0
   ref=np.sum(a[:,0]*a[:,1])/a[:,1].sum();v=C*np.log(d['wavelength_AA']/(ref*(1+self.cfg['z'])))
   use=(v>=self.cfg['window'][0])&(v<=self.cfg['window'][1]);ix=np.flatnonzero(use);v=v[ix]
   dv=float(np.median(np.diff(v)));good=d['valid'][ix]&np.isfinite(d['flux'][ix])&(d['error'][ix]>0)
   row=min(rows,key=lambda r:abs(r['Wavelength']-np.median(d['wavelength_AA'][ix])));nom=C/row['NomResolPower']
   pad=int(np.ceil(6*nom*1.4/dv))+2;grid=v[0]-pad*dv+(np.arange((len(v)+2*pad)*oversample)+.5)*dv/oversample-.5*dv
   self.lines.append(dict(key=k,atom=a,ref=ref,v=v,dv=dv,grid=grid,pad=pad,good=good,flux=d['flux'][ix],error=d['error'][ix],wave=d['wavelength_AA'][ix],source_indices=d['source_indices'][ix],x=(v-v.mean())/100,nominal_fwhm=nom,input_sha256=sha(path)))
   self.slices.append(slice(off,off+good.sum()));off+=good.sum()
  self.ndata=int(off)
 def initial(self,seed=0):
  n=self.ncomp;rng=np.random.default_rng(seed);cc=np.array(self.cfg['centers'])
  if self.target=='J232128-105122':v=np.linspace(4,18,n) if n>1 else np.array([10.])
  else:v=np.linspace(20,53,n) if n>1 else np.array([42.])
  if seed:v+=rng.normal(0,2,n)
  p=list(np.full(n,self.cfg['logtotal']-np.log10(n))+rng.normal(0,.10,n))+list(v)+list(np.log(np.full(n,4.)))
  lo=[10.]*n+[self.cfg['component_bounds'][0]]*n+[np.log(.5)]*n
  hi=[16.]*n+[self.cfg['component_bounds'][1]]*n+[np.log(25.)]*n
  for L in self.lines:
   p += [0,0,0,np.log(L['nominal_fwhm'])]
   lo += [-.05,-.10,-.02,np.log(.6*L['nominal_fwhm'])]
   hi += [.05,.10,.02,np.log(1.4*L['nominal_fwhm'])]
  p += [0]*len(self.shiftkeys);lo += [-1.]*len(self.shiftkeys);hi += [1.]*len(self.shiftkeys)
  return np.array(p),np.array(lo),np.array(hi)
 def evaluate(self,p):
  if self.cache is not None and np.array_equal(p,self.cache[0]):return self.cache[1:]
  n=self.ncomp;N=10**p[:n];vc=p[n:2*n];b=np.exp(p[2*n:3*n]);rr=[];jj=[];profiles=[]
  for ki,L in enumerate(self.lines):
   sh=p[3*n+4*len(self.keys)+self.shiftkeys.index(L['key'])] if L['key'] in self.shiftkeys else 0
   g=L['grid'];tau=np.zeros((len(g),n));tv=np.zeros_like(tau);tb=np.zeros_like(tau)
   for lam,f,gamma,mass,q in L['atom']:
    ratio=np.exp((vc[None,:]+C*np.log(lam/L['ref'])+sh-g[:,None])/C);u=C*(1-ratio)/b[None,:]
    damping=gamma*lam*1e-13/(4*np.pi*b);z=u+1j*damping[None,:];w=wofz(z);wp=-2*z*w+2j/np.sqrt(np.pi)
    amp=TAU*N*f*lam/b;t=w.real*amp[None,:];tau+=t;tv+=amp[None,:]*wp.real*(-ratio/b[None,:]);tb+=-t+amp[None,:]*(-u*wp.real+damping[None,:]*wp.imag)
   intrinsic=np.exp(-tau.sum(axis=1));allarr=np.column_stack([intrinsic,-intrinsic[:,None]*tau*np.log(10),-intrinsic[:,None]*tv,-intrinsic[:,None]*tb,-intrinsic*tv.sum(axis=1)])
   start=3*n+4*ki;c0,c1,zero,lf=p[start:start+4];fw=np.exp(lf);sig=fw/2.354820045/(L['dv']/self.oversample)
   def pix(a,s):return gaussian_filter1d(a,s,axis=0,mode='nearest',truncate=6.).reshape(-1,self.oversample,a.shape[1]).mean(axis=1)[L['pad']:-L['pad']]
   pixarr=pix(allarr,sig);eps=1e-4;dlsf=(pix(intrinsic[:,None],sig*np.exp(eps))[:,0]-pix(intrinsic[:,None],sig*np.exp(-eps))[:,0])/(2*eps)
   cont=1+c0+c1*L['x'];profile=zero+(1-zero)*pixarr[:,0];m=cont*profile;J=np.zeros((len(m),self.npar))
   J[:,:3*n]=(cont*(1-zero))[:,None]*pixarr[:,1:1+3*n];J[:,start]=profile;J[:,start+1]=L['x']*profile;J[:,start+2]=cont*(1-pixarr[:,0]);J[:,start+3]=cont*(1-zero)*dlsf
   if L['key'] in self.shiftkeys:J[:,3*n+4*len(self.keys)+self.shiftkeys.index(L['key'])]=cont*(1-zero)*pixarr[:,-1]
   good=L['good'];rr.append((m[good]-L['flux'][good])/L['error'][good]);jj.append(J[good]/L['error'][good,None]);profiles.append(m)
  self.cache=(p.copy(),np.concatenate(rr),np.vstack(jj),profiles);return self.cache[1:]
 def fun(self,p):return self.evaluate(p)[0]
 def jac(self,p):return self.evaluate(p)[1]
def fit(target,ncomp,free,seed,start=None,suffix='',max_nfev=300,oversample=11):
 m=Model(target,ncomp,free,oversample);p,lo,hi=m.initial(seed)
 if start is not None:
  lookup=dict(zip(start['labels'],start['parameters']));p=np.array([lookup.get(k,v) for k,v in zip(m.labels,p)])
 p=np.clip(p,lo+1e-8,hi-1e-8);initial=p.copy();t=time.time();opt=least_squares(m.fun,p,jac=m.jac,bounds=(lo,hi),x_scale='jac',ftol=2e-8,xtol=1e-8,gtol=1e-6,max_nfev=max_nfev)
 r,J,profiles=m.evaluate(opt.x);chi=float(r@r);ndf=m.ndata-m.npar
 name=f'{target}_n{ncomp}_{"shift" if free else "null"}_s{seed}{suffix}'
 out=dict(name=name,target=target,ncomp=ncomp,free_shifts=free,seed=seed,oversample=oversample,ndata=m.ndata,npar=m.npar,chi2=chi,nominal_ndf=ndf,chi2_per_ndf=chi/ndf,optimizer_success=bool(opt.success),optimizer_status=int(opt.status),message=opt.message,nfev=int(opt.nfev),optimality=float(opt.optimality),seconds=time.time()-t,labels=m.labels,parameters=opt.x.tolist(),initial_parameters=initial.tolist(),bounds_lower=lo.tolist(),bounds_upper=hi.tolist(),active_bounds=[m.labels[i] for i in range(m.npar) if min(opt.x[i]-lo[i],hi[i]-opt.x[i])<1e-5],shifts_m_s={k:1000*float(opt.x[m.labels.index(f'shift_{k}')]) for k in m.shiftkeys},aicc=chi+2*m.npar+2*m.npar*(m.npar+1)/(m.ndata-m.npar-1),source_files={str(DATA/f'{target}_z{m.cfg["z"]:.3f}_FeII{L["key"]}_window.npz'):L['input_sha256'] for L in m.lines},atomic_sha256=sha(ATOMIC),model_script_sha256=sha(__file__),per_line=[dict(line=L['key'],ndata=int(L['good'].sum()),chi2=float(r[sl]@r[sl])) for L,sl in zip(m.lines,m.slices)])
 OUT.mkdir(parents=True,exist_ok=True);(OUT/f'{name}.json').write_text(json.dumps(out,indent=2)+'\n')
 ar=dict(parameters=opt.x,residual=r,jacobian=J)
 for L,prof in zip(m.lines,profiles):
  for key in ['v','wave','flux','error','good','source_indices']:ar[f'{L["key"]}_{key}']=L[key]
  ar[f'{L["key"]}_model']=prof
 np.savez_compressed(OUT/f'{name}.npz',**ar)
 print(name,'chi',chi,'red',chi/ndf,'ok',opt.success,'bound',out['active_bounds'],flush=True)
 return out

def run_target(target):
 records=[]
 for n in [1,2,3,4,6]:
  fits=[fit(target,n,False,s) for s in [0,1,2]]
  chosen=min(fits,key=lambda r:r['chi2'])
  if not chosen['optimizer_success']:chosen=fit(target,n,False,chosen['seed'],chosen,'_continued',600);fits.append(chosen)
  alternatives=[fit(target,n,True,chosen['seed'],chosen,'_fromnull')]
  if not alternatives[0]['optimizer_success']:alternatives.append(fit(target,n,True,chosen['seed'],alternatives[0],'_continued',600))
  alt=min(alternatives,key=lambda r:r['chi2'])
  nullcross=fit(target,n,False,chosen['seed'],alt,'_cross',400);fits.append(nullcross);chosen=min(fits,key=lambda r:r['chi2'])
  records.append(dict(ncomp=n,null=chosen['name'],alternative=alt['name'],chi2_null=chosen['chi2'],chi2_alternative=alt['chi2'],delta_chi2=chosen['chi2']-alt['chi2'],ndata=chosen['ndata'],nominal_ndf_null=chosen['nominal_ndf'],aicc_null=chosen['aicc'],shifts_m_s=alt['shifts_m_s'],both_success=chosen['optimizer_success'] and alt['optimizer_success'],null_active_bounds=chosen['active_bounds'],alternative_active_bounds=alt['active_bounds']))
 (OUT/f'{target}_comparison.json').write_text(json.dumps(dict(configuration=CONFIG[target],comparisons=records,interpretation='Exploratory conventional-profile feasibility. Components/windows selected after inspecting profiles, diagonal pixel errors, unconstrained calibration residuals. No alpha or cosmic-time-law measurement.'),indent=2)+'\n')
 return records
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--target',choices=list(CONFIG));args=a.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 (OUT/'policy.json').write_text(json.dumps(dict(configurations=CONFIG,component_counts=[1,2,3,4,6],start_seeds=[0,1,2],shift_anchor=2374,shift_bounds_km_s=[-1,1],continuum0_bounds=[-.05,.05],continuum1_per100kms_bounds=[-.10,.10],zero_bounds=[-.02,.02],lsf_bounds_relative_to_nominal=[.6,1.4],interpretation='Exploratory fixed windows; no clipping by fitted residual, all source-valid pixels retained. No model selected by an anomalous shift. Fit adequacy evaluated before offsets.'),indent=2)+'\n')
 if args.target:run_target(args.target)
 else:
  with ProcessPoolExecutor(2) as pool:list(pool.map(run_target,CONFIG))
