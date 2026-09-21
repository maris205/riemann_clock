#!/usr/bin/env python3
"""Matched joint native-pixel fits with centered positive asymmetric LSF controls."""
from pathlib import Path
import json,time,hashlib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d
from scipy.linalg import qr,solve_triangular
from scipy.optimize import least_squares
import exposure_analysis as src
ROOT=src.ROOT;OUT=ROOT/'results/order_response';OUT.mkdir(exist_ok=True)
META=json.loads((src.PRO/'metadata.json').read_text());OLD=json.loads((src.OUT/'fit_rows.json').read_text())
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

class JointModel:
 def __init__(self,line=2600,kind='G0',indices=None,frozen=None,step=.025,quadrature=7):
  self.line=line;self.kind=kind;self.asymmetric=kind.startswith('A');self.offset=kind.endswith('1');self.template=src.Template(step=step,quadrature=quadrature);self.frozen=frozen or {};self.indices=list(range(17)) if indices is None else list(indices)
  self.labels=[f'velocity_{i}' for i in self.indices]
  if self.offset and 'order_offset' not in self.frozen:self.labels+=['order_offset']
  self.orders=sorted(set(s['order_index'] for i in self.indices for s in META['exposures'][i]['segments'] if s['line']==line));assert len(self.orders)==2
  for order in self.orders:
   if f'width_{order}' not in self.frozen:self.labels+=[f'width_{order}']
   if self.asymmetric and f'asymmetry_{order}' not in self.frozen:self.labels+=[f'asymmetry_{order}']
  self.blocks=[];self.ndata=0
  for ei in self.indices:
   rec=META['exposures'][ei]
   with np.load(ROOT/rec['array_file']) as arr:
    for s in rec['segments']:
     if s['line']!=line:continue
     pre=s['prefix'];g=arr[pre+'_good'];v=arr[pre+'_v'][g]
     block=dict(exposure=ei,row=s['row'],order=s['order_index'],trace=s['trace_index'],coordinates=self.template.coordinates(arr,pre,line)[g],x=(v-v.mean())/100.,flux=arr[pre+'_flux'][g],error=arr[pre+'_error'][g],velocity=v,pixel=arr[pre+'_pixel'][g]);self.blocks.append(block);self.ndata+=len(v)
  self.npar=len(self.labels);self.nlinear=3*len(self.blocks);self.cache=None
 def initial(self):
  lookup={r['exposure_index']:r['shift_km_s'] for r in OLD if r['configuration']=='primary' and r['line']==self.line}
  p=[];lo=[];hi=[]
  for label in self.labels:
   if label.startswith('velocity_'):p.append(lookup[int(label.split('_')[-1])]);lo.append(-1.5);hi.append(1.5)
   elif label=='order_offset':p.append(0.);lo.append(-.3);hi.append(.3)
   elif label.startswith('width_'):p.append(src.REGIONS[self.line][2]);lo.append(1.6);hi.append(2.6)
   else:p.append(0.);lo.append(-1.);hi.append(1.)
  return np.array(p),np.array(lo),np.array(hi)
 def shape_spline(self,width,asym):
  sigma_eq=width/2.354820045;d=1.2*np.cbrt(asym);sigma=np.sqrt(sigma_eq*sigma_eq-.16*d*d)
  conv=gaussian_filter1d(self.template.intrinsic[self.line],sigma/self.template.step,mode='nearest',truncate=7.)
  sp=CubicSpline(self.template.grid,conv,extrapolate=True)
  values=.8*sp(self.template.grid+.2*d)+.2*sp(self.template.grid-.8*d)
  return CubicSpline(self.template.grid,values,extrapolate=False)
 def evaluate(self,p):
  if self.cache is not None and np.array_equal(p,self.cache[0]):return self.cache[1:]
  pars={**self.frozen,**dict(zip(self.labels,p))};splines={}
  for order in self.orders:
   w=pars[f'width_{order}'];a=pars.get(f'asymmetry_{order}',0.);sp=self.shape_spline(w,a);splines[order]=(sp,sp.derivative())
   for name,h in [('width',1e-4),('asymmetry',1e-4)]:
    label=f'{name}_{order}'
    if label not in self.labels:continue
    plus=self.shape_spline(w+(h if name=='width' else 0),a+(h if name=='asymmetry' else 0));minus=self.shape_spline(w-(h if name=='width' else 0),a-(h if name=='asymmetry' else 0))
    splines[label]=(plus,minus,h)
  residuals=[];jacrows=[];details=[]
  for b in self.blocks:
   order=b['order'];shift=pars[f'velocity_{b["exposure"]}']+(pars.get('order_offset',0.) if order==self.orders[0] else 0.)
   coords=b['coordinates']-shift;sp,dp=splines[order];prof=sp(coords)@self.template.weights/2
   D=np.column_stack([prof,b['x']*prof,np.ones(len(prof))]);X=D/b['error'][:,None];y=b['flux']/b['error'];Q,R=qr(X,mode='economic');beta=solve_triangular(R,Q.T@y);model=D@beta;r=(model-b['flux'])/b['error'];J=np.zeros((len(r),self.npar))
   derivatives={f'velocity_{b["exposure"]}':-(dp(coords)@self.template.weights/2)}
   if self.offset and order==self.orders[0] and 'order_offset' in self.labels:derivatives['order_offset']=derivatives[f'velocity_{b["exposure"]}']
   for name in ['width','asymmetry']:
    label=f'{name}_{order}'
    if label in self.labels:
     plus,minus,h=splines[label];derivatives[label]=((plus(coords)-minus(coords))@self.template.weights/2)/(2*h)
   for label,dprof in derivatives.items():
    dX=np.column_stack([dprof,b['x']*dprof,np.zeros(len(prof))])/b['error'][:,None];z=dX@beta
    # Exact differentiated variable projection: retain residual-dependent term.
    J[:,self.labels.index(label)]=z-Q@(Q.T@z)-Q@solve_triangular(R.T,dX.T@r,lower=True)
   residuals.append(r);jacrows.append(J);details.append(dict(exposure=b['exposure'],row=b['row'],order=order,trace=b['trace'],beta=beta,model=model,residual=r))
  self.cache=(p.copy(),np.concatenate(residuals),np.vstack(jacrows),details);return self.cache[1:]
 def fun(self,p):return self.evaluate(p)[0]
 def jac(self,p):return self.evaluate(p)[1]

def fit(name,line,kind,indices=None,start=None,frozen=None,max_nfev=180):
 m=JointModel(line,kind,indices,frozen);p,lo,hi=m.initial()
 if start:
  p=np.array([start.get(label,value) for label,value in zip(m.labels,p)])
 p=np.clip(p,lo+1e-9,hi-1e-9);t=time.time();opt=least_squares(m.fun,p,jac=m.jac,bounds=(lo,hi),x_scale='jac',ftol=1e-9,xtol=1e-10,gtol=1e-7,max_nfev=max_nfev)
 r,J,details=m.evaluate(opt.x);u,sv,vh=np.linalg.svd(J,full_matrices=False);keep=sv>sv[0]*1e-10;cov=(vh[keep].T/sv[keep]**2)@vh[keep];cor=cov/np.sqrt(np.diag(cov)[:,None]*np.diag(cov)[None,:]);active=[label for label,x,l,h in zip(m.labels,opt.x,lo,hi) if min(x-l,h-x)<1e-5]
 out=dict(name=name,line=line,kind=kind,exposure_indices=m.indices,parameters=dict(zip(m.labels,opt.x.tolist())),frozen=m.frozen,chi2=float(r@r),ndata=m.ndata,npar_nonlinear=m.npar,npar_profiled_linear=m.nlinear,optimizer_success=bool(opt.success),message=opt.message,nfev=int(opt.nfev),optimality=float(opt.optimality),seconds=time.time()-t,labels=m.labels,covariance=cov.tolist(),correlation=cor.tolist(),singular_values=sv.tolist(),condition_number=float(sv[0]/sv[-1]),active_bounds=active,start_parameters=p.tolist(),per_row=[dict(exposure=d['exposure'],row=d['row'],order=d['order'],trace=d['trace'],coefficients=d['beta'].tolist(),chi2=float(d['residual']@d['residual']),npix=len(d['residual'])) for d in details])
 if 'order_offset' in m.labels:
  j=m.labels.index('order_offset');out['order_offset_m_s']=float(opt.x[j]*1000);out['order_offset_sigma_m_s']=float(np.sqrt(cov[j,j])*1000);out['offset_correlations']={l:float(cor[j,i]) for i,l in enumerate(m.labels) if l.startswith(('width_','asymmetry_'))}
 arrays=dict(parameters=opt.x,residual=r,jacobian=J)
 for d in details:
  pre=f'exp{d["exposure"]}_row{d["row"]}';arrays[pre+'_model']=d['model'];arrays[pre+'_residual']=d['residual']
 (OUT/f'{name}.json').write_text(json.dumps(out,indent=2)+'\n');np.savez_compressed(OUT/f'{name}.npz',**arrays)
 print(name,'chi2',out['chi2'],'offset',out.get('order_offset_m_s'),'success',opt.success,'active',active,'seconds',round(out['seconds'],2),flush=True)
 return out

def run_block(line,indices,label):
 selected={}
 for kind in ['G0','G1','A0','A1']:
  candidates=[]
  starts=[(0.,0.)] if kind.startswith('G') else [(0.,0.),(.5,-.5),(-.5,.5)]
  for si,seeds in enumerate(starts):
   start={}
   if kind=='G1':start.update(selected['G0']['parameters'])
   if kind.startswith('A'):start.update(selected['G'+kind[-1]]['parameters'])
   orders=JointModel(line,kind,indices).orders
   if kind.startswith('A'):
    for order,value in zip(orders,seeds):start[f'asymmetry_{order}']=value
   candidates.append(fit(f'{label}_{line}_{kind}_start{si}',line,kind,indices,start))
  successful=[r for r in candidates if r['optimizer_success']]
  selected[kind]=min(successful or candidates,key=lambda r:r['chi2'])
  (OUT/f'{label}_{line}_selected.json').write_text(json.dumps(selected,indent=2)+'\n')
 return selected

def main():
 allout={}
 for line in [2600,2382]:
  full=run_block(line,list(range(17)),'full');train=run_block(line,list(range(9)),'train2018');heldout={}
  for kind,source in train.items():
   frozen={k:v for k,v in source['parameters'].items() if not k.startswith('velocity_')}
   heldout[kind]=fit(f'heldout_{line}_{kind}',line,kind,list(range(9,17)),frozen=frozen)
  out=dict(full=full,train2018=train,heldout2019_2020=heldout)
  allout[str(line)]=out;(OUT/'model_summary.json').write_text(json.dumps(allout,indent=2)+'\n')
 (OUT/'model_run_manifest.json').write_text(json.dumps(dict(code_sha256=SHA(__file__),protocol_sha256=SHA(OUT/'model_protocol.json'),source_metadata_sha256=SHA(src.PRO/'metadata.json'),source_template_sha256=SHA(src.TEMPLATE_PATH),source_previous_fits_sha256=SHA(src.OUT/'fit_rows.json'),completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())),indent=2)+'\n')
if __name__=='__main__':main()
