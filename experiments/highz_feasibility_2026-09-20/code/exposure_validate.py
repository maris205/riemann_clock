#!/usr/bin/env python3
"""Independent reconstruction and numerical audit of the native exposure diagnostic.
Uses raw FITS + original UPL rather than the prepared mask artifact for extraction checks;
QR variable projection, centered derivative steps and independent GLS summary equations.
"""
from pathlib import Path
import sys,json,re,hashlib,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'.deps'))
import numpy as np
from astropy.io import fits
from scipy.linalg import qr,solve_triangular
from scipy.stats import chi2
from scipy.special import voigt_profile
from scipy.signal import fftconvolve
import exposure_analysis as a
OUT=ROOT/'results/exposures';checks=[];metrics={}
def check(name,yes,details=None):
 checks.append(dict(name=name,pass_=bool(yes),details=details))
def close(name,x,y,rtol=1e-8,atol=1e-10):
 xx=np.asarray(x);yy=np.asarray(y);err=float(np.max(abs(xx-yy))) if xx.size else 0.
 check(name,np.allclose(xx,yy,rtol=rtol,atol=atol),{'maximum_absolute_error':err})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
meta=json.loads((a.PRO/'metadata.json').read_text());records=meta['exposures'];rows=json.loads((OUT/'fit_rows.json').read_text());summary=json.loads((OUT/'summary.json').read_text());config=json.loads((OUT/'configuration.json').read_text())
check('analysis source matches saved configuration',config['script_sha256']==sha(a.__file__))
check('template source matches saved configuration',config['template_sha256']==sha(a.TEMPLATE_PATH))
manifest=json.loads((a.RAW/'manifest.json').read_text());upl=(ROOT/'data/raw/espresso_null/hes0515m4414.upl').read_text();filemap=dict((int(m[0])-1,m[1]) for m in re.findall(r'(?m)^\s*(\d+)_FLUX = (.+)$',upl));clips=[]
for b in re.split(r'(?m)^(?=\d+_ACTN)',upl)[1:]:
 typ=int(re.search(r'_ACTN = (\d+)',b)[1])
 if typ!=1:continue
 co=[float(x) for x in re.search(r'_CORD =\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)',b).groups()]
 se,ro=[int(x)-1 for x in re.search(r'_SPEC = (\d+)\s+\d+_ORDR = (\d+)',b).groups()];clips.append((se,ro,*co))
check('UPL has no unclip actions',not re.search(r'_ACTN = [34]\s',upl))
check('17 raw exposures',len(records)==len(manifest['files'])==17)
raw_audits=[];arrays={};maxwidth=0.;allflags={};counts=0;atmoscounts=0
atmosphere=np.loadtxt(ROOT/'data/raw/espresso_null/atmomask.dat')
check('atmospheric source hash',meta['atmosphere_mask_sha256']==sha(ROOT/'data/raw/espresso_null/atmomask.dat'))
for rec,it in zip(records,manifest['files']):
 ei=rec['index'];p=ROOT/rec['source_file'];data=p.read_bytes();check(f'raw{ei} SHA256',hashlib.sha256(data).hexdigest()==it['sha256']);check(f'raw{ei} git blob',hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()==it['git_blob_sha1']);del data
 check(f'raw{ei} UPL exposure mapping',filemap[ei]==p.name)
 arr=dict(np.load(ROOT/rec['array_file']));arrays[ei]=arr
 with fits.open(p,memmap=False) as h:
  w=h['WAVEDATA_VAC_BARY'].data;fl=h['SCIDATA'].data;er=h['ERRDATA'].data;quality=h['QUALDATA'].data;dll=h['DLLDATA_VAC_BARY'].data
  check(f'raw{ei} schema',w.shape==fl.shape==er.shape==quality.shape==dll.shape==(170,9111))
  check(f'raw{ei} DRS version',h[0].header['ESO PRO REC1 PIPE ID']=='espdr/2.2.3')
  check(f'raw{ei} sky product',h[0].header['ESO PRO CATG']=='S2D_SKYSUB_A')
  for val,n in zip(*np.unique(quality,return_counts=True)):allflags[str(int(val))]=allflags.get(str(int(val)),0)+int(n)
  for s in rec['segments']:
   pre=s['prefix'];row=s['row'];k=s['line'];lo,hi,_=a.REGIONS[k];idx=np.flatnonzero((w[row]>=lo)&(w[row]<=hi));center=w[row,idx];left=.5*(w[row,idx-1]+center);right=.5*(center+w[row,idx+1])
   close(f'{ei}/{pre} pixels',arr[pre+'_pixel'],idx);close(f'{ei}/{pre} WAVE',arr[pre+'_wave'],center);close(f'{ei}/{pre} left',arr[pre+'_left'],left);close(f'{ei}/{pre} right',arr[pre+'_right'],right)
   # Independent nearest neighbor lookup through full distances, not searchsorted.
   minc,maxc=np.searchsorted(a.COADD_WAVE,[lo-.1,hi+.1]);cw=a.COADD_WAVE[minc:maxc];near=minc+np.argmin(abs(cw[:,None]-center),axis=0)
   good=(quality[row,idx]==0)&np.isfinite(fl[row,idx])&np.isfinite(er[row,idx])&(er[row,idx]>0)&(a.COADD_DATA[4,near]==1);before=good.copy()
   for exposure,trace,l,r,ylo,yhi in clips:
    if exposure==ei and trace==row and l<=hi and r>=lo:
     check(f'{ei}/{pre} clip flux range irrelevant {l}',np.min(fl[row,idx])>ylo and np.max(fl[row,idx])<yhi)
     good &= ~((left<=r)&(right>=l))
   afterupl=good.copy()
   berv=float(h[0].header['ESO QC BERV']);pad=meta['atmosphere_padding_km_s']
   shifted=atmosphere*(1+(berv+np.array([-pad,pad]))/a.C)
   atmospheric_overlap=((right[:,None]>=shifted[:,0])&(left[:,None]<=shifted[:,1])).any(axis=1)
   good &= ~atmospheric_overlap
   check(f'{ei}/{pre} exact fixed mask',np.array_equal(good,arr[pre+'_good']))
   check(f'{ei}/{pre} mask counts',int(good.sum())==s['ngood'] and int((before&~afterupl).sum())==s['additional_upl_masked_pixels'] and int((afterupl&~good).sum())==s['additional_atmospheric_masked_pixels'])
   scale=float(np.quantile(fl[row,idx][good],.9));close(f'{ei}/{pre} scale',s['scale_counts'],scale)
   close(f'{ei}/{pre} flux',arr[pre+'_flux'],fl[row,idx].astype(float)/scale);close(f'{ei}/{pre} sigma',arr[pre+'_error'],er[row,idx].astype(float)/scale)
   ref=np.average(a.ATOMS[k][:,0],weights=a.ATOMS[k][:,1]);close(f'{ei}/{pre} velocity',arr[pre+'_v'],a.C*np.log(center/(ref*(1+a.Z))))
   check(f'{ei}/{pre} row trace mapping',s['order_index']==row//2 and s['trace_index']==row%2)
   maxwidth=max(maxwidth,float(np.max(abs((right-left)/dll[row,idx]-1))));counts+=int((before&~afterupl).sum());atmoscounts+=int((afterupl&~good).sum())
metrics.update(raw_quality_counts=allflags,additional_manual_mask_pixels=counts,additional_atmosphere_mask_pixels=atmoscounts,max_midpoint_vs_DLL_relative_width_difference=maxwidth)
T={k:a.Template(k) for k in ['null_cross','published']}

def replay(fit,p,freeze=None,ret=False):
 rec=records[fit['exposure_index']];arr=arrays[rec['index']];t=T[fit['template']];k=fit['line'];fwhm=p[1] if fit['free_lsf'] else a.REGIONS[k][2];sp=t.spline(k,fwhm);res=[];dets=[];nuis=[]
 for s in a.choose_segments(rec,k,fit['selection']):
  pre=s['prefix'];g=arr[pre+'_good'];vel=arr[pre+'_v'][g];x=(vel-vel.mean())/100
  v=t.coordinates(arr,pre,k)[g];prof=sp(v-p[0])@t.weights/2;D=np.column_stack([prof,x*prof,np.ones(len(x))]);err=arr[pre+'_error'][g];y=arr[pre+'_flux'][g];X=D/err[:,None];Q,R=qr(X,mode='economic')
  beta=solve_triangular(R,Q.T@(y/err)) if freeze is None else freeze[len(dets)]
  model=D@beta;rr=(model-y)/err;res.append(rr);dets.append((s,beta,model,rr,Q));nuis.append(Q)
 return (np.concatenate(res),dets) if ret else np.concatenate(res)
maxcov=0.;maxproj=0.;maxgrad=0.;bounds=[]
for i,r in enumerate(rows):
 p=np.array(r['parameters']);rr,dets=replay(r,p,ret=True);name=f"{r['configuration']}/{r['exposure_index']}/{r['line']}"
 close(name+' chi2',r['chi2'],rr@rr,rtol=1e-10,atol=1e-7)
 check(name+' N and P',len(rr)==r['ndata'] and len(p)+3*len(dets)==r['npar'])
 for (_,b,m,res,Q),old in zip(dets,r['rows']):
  close(name+f"/{old['row']} model",old['model'],m,atol=1e-9);close(name+f"/{old['row']} residual",old['residual'],res,atol=1e-8)
  close(name+f"/{old['row']} nuisance",[old['continuum_amplitude'],old['continuum_slope'],old['additive_zero']],b,atol=1e-9)
 check(name+' optimizer success',r['optimizer_success'])
 # Centered derivative of independently profiled residual at two steps.
 def derivative(step,freeze=None):
  cols=[]
  for j in range(len(p)):
   d=np.zeros(len(p));d[j]=step;cols.append((replay(r,p+d,freeze=freeze)-replay(r,p-d,freeze=freeze))/(2*step))
  return np.column_stack(cols)
 J=derivative(1e-5);J2=derivative(3e-5);close(name+' derivative step stability',J,J2,rtol=1e-4,atol=1e-6)
 cov=np.linalg.inv(J.T@J);oldcov=np.array(r['covariance']);rel=float(np.max(abs(cov-oldcov))/np.max(abs(cov)));maxcov=max(maxcov,rel)
 check(name+' covariance replay',rel<.003,{'relative_error':rel})
 # Different nuisance-projected Fisher construction with beta held fixed.
 JJ=derivative(1e-5,freeze=[d[1] for d in dets]);off=0
 for _,_,_,res,Q in dets:
  sl=slice(off,off+len(res));JJ[sl]-=Q@(Q.T@JJ[sl]);off+=len(res)
 cvproj=np.linalg.inv(JJ.T@JJ);maxproj=max(maxproj,abs(float(np.sqrt(cvproj[0,0]/cov[0,0]))-1))
 # Gradient metric is diagnostic; nominal optimizer success is not global uniqueness.
 grad=float(np.max(abs(J.T@rr)));maxgrad=max(maxgrad,grad)
 if r['bound_active']:bounds.append(name)
metrics.update(max_relative_covariance_replay_error=maxcov,max_profiled_vs_frozen_nuisance_projected_sigma_fraction=maxproj,max_absolute_profile_gradient=maxgrad,bound_active_fits=bounds)
for name,s in summary.items():
 points=s['points'];fits=[r for r in rows if r['configuration']==name]
 for p in points:
  fr={r['line']:r for r in fits if r['exposure_index']==p['exposure_index']};v=np.array([fr[k]['shift_km_s'] for k in a.KEYS]);var=np.array([fr[k]['shift_sigma_km_s']**2 for k in a.KEYS]);A=np.array([[-1,1,0],[-1,0,1]])
  close(f'{name}/{p["exposure_index"]} relative',p['relative_m_s'],A@v*1000)
  close(f'{name}/{p["exposure_index"]} anchor covariance',p['covariance_m2_s2'],1e6*A@np.diag(var)@A.T)
 for part in ['all_exposures','chronological_2018','chronological_2019_2020']:
  out=s[part];ps=[p for p in points if p['exposure_index'] in out['exposure_indices']];n=len(ps);B=np.tile(np.eye(2),(n,1));V=np.zeros((2*n,2*n));y=np.concatenate([p['relative_m_s'] for p in ps])
  for j,p in enumerate(ps):V[2*j:2*j+2,2*j:2*j+2]=p['covariance_m2_s2']
  Vinv=np.linalg.inv(V);cov=np.linalg.inv(B.T@Vinv@B);m=cov@B.T@Vinv@y;d=y-B@m;q=d@Vinv@d
  close(f'{name}/{part} mean',out['mean_relative_m_s'],m);close(f'{name}/{part} covariance',out['conditional_covariance_m2_s2'],cov);close(f'{name}/{part} heterogeneity',out['heterogeneity_chi2'],q)
 early=set(s['chronological_2018']['exposure_indices']);late=set(s['chronological_2019_2020']['exposure_indices']);check(name+' groups disjoint complete',len(early)==9 and len(late)==8 and not early&late and early|late==set(range(17)))
 d=np.array(s['chronological_2018']['mean_relative_m_s'])-s['chronological_2019_2020']['mean_relative_m_s'];cv=np.array(s['chronological_2018']['conditional_covariance_m2_s2'])+s['chronological_2019_2020']['conditional_covariance_m2_s2'];q=d@np.linalg.solve(cv,d)
 close(name+' chronology chi2',s['chronological_difference']['chi2'],q);close(name+' chronology p',s['chronological_difference']['conditional_p'],chi2.sf(q,2))
# Intrinsic profile independently computed using scipy.voigt_profile and explicit sum.
for kind,t in T.items():
 for k in a.KEYS:
  v=t.grid[::137];atom=a.ATOMS[k];ref=np.average(atom[:,0],weights=atom[:,1]);tau=np.zeros(len(v))
  for logN,vc,b in t.comp:
   for lam,f,gamma,_,_ in atom:
    u=a.C*(1-np.exp((vc+a.C*np.log(lam/ref)-v)/a.C));gammavel=gamma*lam*1e-13/(4*np.pi)
    # wofz(u/b+i*g/b).real = sqrt(pi)*b*Voigt(u,b/sqrt(2),g).
    tau+=a.core.TAU_CONSTANT*10**logN*f*lam*np.sqrt(np.pi)*voigt_profile(u,b/np.sqrt(2),gammavel)
  close(f'{kind}/{k} independent Voigt normalization',t.intrinsic[k][::137],np.exp(-tau),rtol=1e-12,atol=1e-12)
# Refine template grid 2x and quadrature 3x for representative primary, free-LSF,
# and original published-template fits. Selection is fixed to exposure 0.
refined={k:a.Template(k,step=.0125,quadrature=21) for k in ['null_cross','published']};refine=[]
for r in rows:
 if r['exposure_index']!=0 or r['configuration'] not in ['primary','free_lsf','published_template']:continue
 rec=records[0];new=a.fit_line(rec,arrays[0],r['line'],refined[r['template']],r['selection'],r['free_lsf']);delta=1000*(new['shift_km_s']-r['shift_km_s']);dc=new['chi2']-r['chi2'];refine.append(dict(configuration=r['configuration'],line=r['line'],delta_shift_m_s=delta,delta_chi2=dc));check(f"refine {r['configuration']}/{r['line']}",abs(delta)<.1 and abs(dc)<.001,refine[-1])
metrics['refinement']=refine
# Independent paired-difference checks. Orders share exactly the same 2374 anchor,
# while target orders and trace subsets must comprise distinct native rows.
diagnostics=json.loads((OUT/'diagnostics.json').read_text())
for field,first,second,cancel in [('trace_difference','trace0','trace1',False),('order_difference','lower_order','upper_order',True),('order_difference_free_lsf','lower_order_free_lsf','upper_order_free_lsf',True)]:
 ps=[]
 for ei in range(17):
  fr={c:{r['line']:r for r in rows if r['configuration']==c and r['exposure_index']==ei} for c in [first,second]}
  if cancel:
   close(f'{field}/{ei} same anchor shift',fr[first][2374]['shift_km_s'],fr[second][2374]['shift_km_s'],atol=0,rtol=0)
   check(f'{field}/{ei} same anchor data',[d['row'] for d in fr[first][2374]['rows']]==[d['row'] for d in fr[second][2374]['rows']])
  for k in ([2382,2600] if cancel else a.KEYS):
   one={d['row'] for d in fr[first][k]['rows']};two={d['row'] for d in fr[second][k]['rows']};check(f'{field}/{ei}/{k} native subsets disjoint',not one&two)
  v=np.array([[fr[c][k]['shift_km_s'] for k in a.KEYS] for c in [first,second]]);vv=np.array([[fr[c][k]['shift_sigma_km_s']**2 for k in a.KEYS] for c in [first,second]])
  y=(v[0,1:]-v[1,1:])*1000;cv=np.diag(vv[:,1:].sum(axis=0))*1e6
  if not cancel:y-=1000*(v[0,0]-v[1,0]);cv+=1e6*vv[:,0].sum()
  old=diagnostics[field]['points'][ei];close(f'{field}/{ei} contrast',old['relative_m_s'],y);close(f'{field}/{ei} covariance',old['covariance_m2_s2'],cv);ps.append((y,cv))
 W=sum(np.linalg.inv(cv) for _,cv in ps);cv=np.linalg.inv(W);mean=cv@sum(np.linalg.solve(v,y) for y,v in ps);q=mean@np.linalg.solve(cv,mean)
 close(field+' weighted contrast',diagnostics[field]['mean_relative_m_s'],mean);close(field+' weighted covariance',diagnostics[field]['conditional_covariance_m2_s2'],cv);close(field+' joint zero test',diagnostics[field]['zero_difference_chi2'],q)
result=dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),pass_=all(x['pass_'] for x in checks),n_checks=len(checks),n_pass=sum(x['pass_'] for x in checks),n_fail=sum(not x['pass_'] for x in checks),metrics=metrics,scope='Implementation/reproducibility audit; not verification of conditional error model, global fit uniqueness, shared-template independence or new physics.',hashes={p:sha(ROOT/p) for p in ['code/exposure_analysis.py','data/processed/exposures/metadata.json','results/exposures/configuration.json','results/exposures/fit_rows.json','results/exposures/summary.json','results/exposures/diagnostics.json']},checks=checks)
(OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['pass_','n_checks','n_pass','n_fail']},indent=2));print(json.dumps(metrics,indent=2));print('Failures:',[c for c in checks if not c['pass_']])
