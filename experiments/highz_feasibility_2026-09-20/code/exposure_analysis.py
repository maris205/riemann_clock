#!/usr/bin/env python3
"""Native-pixel fixed-Voigt-template consistency checks in 17 ESPRESSO exposures.
No cosmic-time inference and no independent remeasurement of the gas model.
"""
from pathlib import Path
import sys,json,hashlib,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'.deps'))
import numpy as np
from astropy.io import fits
from scipy.special import wofz
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares
from scipy.stats import chi2
import espresso_conventional_null as core
RAW=ROOT/'data/raw/exposures';PRO=ROOT/'data/processed/exposures';OUT=ROOT/'results/exposures'
for d in [PRO,OUT]:d.mkdir(exist_ok=True,parents=True)
KEYS=[2374,2382,2600];C=core.C;Z=core.ZREF
COMP,ATOMS,REGIONS,COADD_WAVE,COADD_DATA,_=core.load_sources()
TEMPLATE_PATH=ROOT/'results/espresso_null/null_cross.json'
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def extract():
 records=[];manifest=json.loads((RAW/'manifest.json').read_text())
 masks=json.loads((RAW/'provenance_upl_line_mask_actions.json').read_text())['actions']
 atmosphere=np.loadtxt(ROOT/'data/raw/espresso_null/atmomask.dat');atmosphere_padding_km_s=.25
 for ei,item in enumerate(manifest['files']):
  path=ROOT/item['local_file']; arr={};segments=[]
  with fits.open(path,memmap=False) as h:
   hdr=h[0].header;w=h['WAVEDATA_VAC_BARY'].data;f=h['SCIDATA'].data;e=h['ERRDATA'].data;q=h['QUALDATA'].data
   for k in KEYS:
    ref=np.average(ATOMS[k][:,0],weights=ATOMS[k][:,1]);lo,hi,_=REGIONS[k]
    for row in range(w.shape[0]):
     ix=np.flatnonzero((w[row]>=lo)&(w[row]<=hi))
     if len(ix)<50:continue
     good=(q[row,ix]==0)&np.isfinite(f[row,ix])&np.isfinite(e[row,ix])&(e[row,ix]>0)
     # Transfer only coadd's broad published region mask, not individual pixel clipping.
     nearest=np.searchsorted(COADD_WAVE,w[row,ix]);nearest=np.clip(nearest,1,len(COADD_WAVE)-1)
     nearest-=abs(COADD_WAVE[nearest-1]-w[row,ix])<abs(COADD_WAVE[nearest]-w[row,ix])
     good &= COADD_DATA[4,nearest]==1
     if good.sum()<50:continue
     center=w[row,ix];left=(w[row,ix-1]+center)/2;right=(center+w[row,ix+1])/2
     good_before_upl=good.copy()
     for mask in masks:
      if mask['exposure_one_based']==ei+1 and mask['row_one_based']==row+1:
       good &= ~((right>=mask['wavelength_lower_A'])&(left<=mask['wavelength_upper_A']))
     good_after_upl=good.copy()
     for atmlo,atmhi in atmosphere:
      barylo=atmlo*(1+(float(hdr['ESO QC BERV'])-atmosphere_padding_km_s)/C)
      baryhi=atmhi*(1+(float(hdr['ESO QC BERV'])+atmosphere_padding_km_s)/C)
      good &= ~((right>=barylo)&(left<=baryhi))
     if good.sum()<50:continue
     prefix=f'{k}_row{row}';scale=float(np.quantile(f[row,ix][good],.9))
     if scale<=0:continue
     arr[prefix+'_wave']=center;arr[prefix+'_left']=left;arr[prefix+'_right']=right
     arr[prefix+'_v']=C*np.log(center/(ref*(1+Z)))
     arr[prefix+'_flux']=f[row,ix].astype(float)/scale;arr[prefix+'_error']=e[row,ix].astype(float)/scale
     arr[prefix+'_quality']=q[row,ix];arr[prefix+'_good']=good;arr[prefix+'_pixel']=ix
     segments.append(dict(prefix=prefix,line=k,row=row,order_index=row//2,trace_index=row%2,npix=len(ix),ngood=int(good.sum()),additional_upl_masked_pixels=int((good_before_upl&~good_after_upl).sum()),additional_atmospheric_masked_pixels=int((good_after_upl&~good).sum()),scale_counts=scale))
   npz=PRO/f'exposure_{ei:02d}.npz';np.savez_compressed(npz,**arr)
   rec=dict(index=ei,date_obs=hdr['DATE-OBS'],mjd_obs=hdr['MJD-OBS'],exptime_s=hdr['EXPTIME'],berv_km_s=hdr.get('ESO QC BERV'),source_file=item['local_file'],source_sha256=item['sha256'],array_file=str(npz.relative_to(ROOT)),segments=segments)
   records.append(rec)
 metadata=dict(atmosphere_mask_sha256=SHA(ROOT/'data/raw/espresso_null/atmomask.dat'),atmosphere_padding_km_s=atmosphere_padding_km_s,upl_masks_sha256=SHA(RAW/'provenance_upl_line_mask_actions.json'),exposures=records,lines=KEYS,source_manifest_sha256=SHA(RAW/'manifest.json'),native_pixels=True,wavelength='WAVEDATA_VAC_BARY: vacuum Angstrom already barycentric; no second BERV correction',error='ERRDATA is one-sigma; no rescaling by fitted chi-square',mask='quality==0, finite flux and sigma>0, exact published wavelength regions, nearest coadd status==1; published UPL exposure/order clipping rectangles conservatively transferred by native-bin overlap (not exact redispersed UPL replay); atmospheric intervals shifted by pipeline BERV with conservative +/-0.25km/s margin and transferred by native-bin overlap',pixel_integration='midpoints of adjacent native wavelength samples, integration in wavelength')
 (PRO/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n');return metadata

class Template:
 def __init__(self,kind='null_cross',step=.025,quadrature=7):
  self.kind=kind;self.step=step;self.quadrature=quadrature
  comp=COMP.copy()
  if kind=='null_cross':
   old=json.loads(TEMPLATE_PATH.read_text());p=np.array(old['parameters']);comp=np.column_stack([p[:45],p[45:90],np.exp(p[90:135])])
  self.comp=comp;self.grid=np.arange(-110,190+step/2,step);self.cache={};self.intrinsic={}
  for k in KEYS:
   atom=ATOMS[k];ref=np.average(atom[:,0],weights=atom[:,1]);tau=np.zeros(len(self.grid))
   for lam,f,gamma,mass,q in atom:
    b=comp[:,2];ratio=np.exp((comp[:,1][None,:]+C*np.log(lam/ref)-self.grid[:,None])/C)
    u=C*(1-ratio)/b[None,:];a=gamma*lam*1e-13/(4*np.pi*b)
    tau+=(wofz(u+1j*a).real*(core.TAU_CONSTANT*10**comp[:,0]*f*lam/b)[None,:]).sum(axis=1)
   self.intrinsic[k]=np.exp(-tau)
  self.nodes,self.weights=np.polynomial.legendre.leggauss(quadrature)
 def spline(self,k,fwhm):
  key=(k,float(fwhm))
  if key not in self.cache:
   conv=gaussian_filter1d(self.intrinsic[k],fwhm/2.354820045/self.step,mode='nearest',truncate=6)
   if len(self.cache)>100:self.cache.clear()
   self.cache[key]=CubicSpline(self.grid,conv,extrapolate=False)
  return self.cache[key]
 def coordinates(self,arr,prefix,k):
  ref=np.average(ATOMS[k][:,0],weights=ATOMS[k][:,1]);left=arr[prefix+'_left'];right=arr[prefix+'_right']
  wave=(left[:,None]+right[:,None])/2+(right-left)[:,None]*self.nodes/2
  return C*np.log(wave/(ref*(1+Z)))


def choose_segments(record,k,selection):
 seg=[s for s in record['segments'] if s['line']==k]
 if selection.startswith('trace'):seg=[s for s in seg if s['trace_index']==int(selection[-1])]
 if selection in ['lower_order','upper_order']:
  order=min(s['order_index'] for s in seg) if selection=='lower_order' else max(s['order_index'] for s in seg)
  seg=[s for s in seg if s['order_index']==order]
 return seg

def fit_line(record,arr,k,template,selection='all',free_lsf=False):
 seg=choose_segments(record,k,selection);blocks=[];ndata=0
 for s in seg:
  pre=s['prefix'];good=arr[pre+'_good'];v=arr[pre+'_v'][good];flux=arr[pre+'_flux'][good];err=arr[pre+'_error'][good]
  blocks.append(dict(meta=s,velocity_nodes=template.coordinates(arr,pre,k)[good],x=(v-v.mean())/100,flux=flux,error=err));ndata+=len(v)
 nominal=REGIONS[k][2]
 def compute(p,return_details=False):
  shift=p[0];fwhm=p[1] if free_lsf else nominal;sp=template.spline(k,fwhm);rr=[];details=[]
  for b in blocks:
   prof=sp(b['velocity_nodes']-shift)@template.weights/2
   # Independent multiplicative continuum amplitude/slope and additive zero for each row.
   design=np.column_stack([prof,b['x']*prof,np.ones(len(prof))]);X=design/b['error'][:,None];y=b['flux']/b['error']
   beta=np.linalg.lstsq(X,y,rcond=1e-12)[0];model=design@beta;r=(model-b['flux'])/b['error'];rr.append(r)
   if return_details:details.append(dict(row=b['meta']['row'],continuum_amplitude=float(beta[0]),continuum_slope=float(beta[1]),additive_zero=float(beta[2]),chi2=float(r@r),ndata=len(r),model=model.tolist(),residual=r.tolist()))
  return (np.concatenate(rr),details) if return_details else np.concatenate(rr)
 p0=[0.,nominal] if free_lsf else [0.];low=[-1.5,1.6] if free_lsf else [-1.5];high=[1.5,2.6] if free_lsf else [1.5]
 opt=least_squares(compute,p0,bounds=(low,high),diff_step=1e-4,ftol=1e-10,xtol=1e-10,gtol=1e-8,max_nfev=100)
 residual,details=compute(opt.x,True);cov=np.linalg.inv(opt.jac.T@opt.jac)
 out=dict(exposure_index=record['index'],date_obs=record['date_obs'],line=k,selection=selection,template=template.kind,free_lsf=free_lsf,shift_km_s=float(opt.x[0]),shift_sigma_km_s=float(np.sqrt(cov[0,0])),fwhm_km_s=float(opt.x[1] if free_lsf else nominal),covariance=cov.tolist(),chi2=float(residual@residual),ndata=ndata,npar=len(opt.x)+3*len(seg),optimizer_success=bool(opt.success),optimality=float(opt.optimality),nfev=int(opt.nfev),bound_active=bool(np.any((opt.x-np.array(low)<1e-5)|(np.array(high)-opt.x<1e-5))),parameters=opt.x.tolist(),rows=details)
 return out

def relative(rows):
 by={r['line']:r for r in rows};v=np.array([by[k]['shift_km_s'] for k in KEYS]);var=np.array([by[k]['shift_sigma_km_s']**2 for k in KEYS]);A=np.array([[-1,1,0],[-1,0,1]])
 return dict(exposure_index=rows[0]['exposure_index'],date_obs=rows[0]['date_obs'],common_velocity_km_s=float(v[0]),relative_m_s=(1000*A@v).tolist(),covariance_m2_s2=(1e6*A@np.diag(var)@A.T).tolist())

def combine(points):
 ys=np.array([p['relative_m_s'] for p in points]);W=np.array([np.linalg.inv(p['covariance_m2_s2']) for p in points]);cov=np.linalg.inv(W.sum(axis=0));mean=cov@np.einsum('nij,nj->i',W,ys);d=ys-mean;q=float(np.einsum('ni,nij,nj',d,W,d));df=2*len(points)-2
 return dict(n_exposures=len(points),exposure_indices=[p['exposure_index'] for p in points],mean_relative_m_s=mean.tolist(),conditional_covariance_m2_s2=cov.tolist(),conditional_sigma_m_s=np.sqrt(np.diag(cov)).tolist(),heterogeneity_chi2=q,heterogeneity_ndf=df,conditional_heterogeneity_p=float(chi2.sf(q,df)) if df else None)

def run():
 meta=extract();template=Template();published=Template('published');allrows=[];collections={}
 configs=[('primary','all',False,template),('free_lsf','all',True,template),('trace0','trace0',False,template),('trace1','trace1',False,template),('lower_order','lower_order',False,template),('upper_order','upper_order',False,template),('published_template','all',False,published),('lower_order_free_lsf','lower_order',True,template),('upper_order_free_lsf','upper_order',True,template)]
 for name,sel,lsf,t in configs:
  points=[]
  for record in meta['exposures']:
   with np.load(ROOT/record['array_file']) as arr:
    rows=[fit_line(record,arr,k,t,sel,lsf) for k in KEYS]
   allrows.extend([dict(configuration=name,**r) for r in rows]);points.append(relative(rows))
   print(name,record['index'],np.round(points[-1]['relative_m_s'],1),flush=True)
  early=[p for p in points if p['date_obs'].startswith('2018')];late=[p for p in points if not p['date_obs'].startswith('2018')]
  a,b=combine(early),combine(late);d=np.array(a['mean_relative_m_s'])-b['mean_relative_m_s'];cov=np.array(a['conditional_covariance_m2_s2'])+b['conditional_covariance_m2_s2'];q=float(d@np.linalg.solve(cov,d))
  collections[name]=dict(points=points,all_exposures=combine(points),chronological_2018=a,chronological_2019_2020=b,chronological_difference=dict(early_minus_late_m_s=d.tolist(),conditional_covariance_m2_s2=cov.tolist(),chi2=q,ndf=2,conditional_p=float(chi2.sf(q,2))))
  (OUT/'fit_rows.json').write_text(json.dumps(allrows,indent=2)+'\n');(OUT/'summary.json').write_text(json.dumps(collections,indent=2)+'\n')
 config=dict(completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),template_source=str(TEMPLATE_PATH.relative_to(ROOT)),template_sha256=SHA(TEMPLATE_PATH),script_sha256=SHA(__file__),quadrature=template.quadrature,grid_step_km_s=template.step,selected_lines=KEYS,configurations=[c[0] for c in configs],interpretation='Fixed shared gas-template native-pixel consistency diagnostic, not full gas refit or independent astrophysical confirmation; template estimated from same coadd. Conditional diagonal ERRDATA covariance, no empirical extracted-pixel covariance available. Exposure/epoch differences test observing systematics, not cosmic time law.')
 (OUT/'configuration.json').write_text(json.dumps(config,indent=2)+'\n')
if __name__=='__main__':run()
