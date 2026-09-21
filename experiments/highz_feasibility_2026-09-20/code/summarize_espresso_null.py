#!/usr/bin/env python3
"""Build transparent comparison tables and actual-data figures from saved fits."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/espresso_null';FIG=ROOT/'figures'
PAIRS=[('primary','null','alternative'),('expected_fluctuation','expected_fluctuation_null','expected_fluctuation_alternative'),('ar1_rho03','ar1_rho03_null','ar1_rho03_alternative'),('lsf_minus3pct','lsf_minus3pct_null','lsf_minus3pct_alternative'),('lsf_plus3pct','lsf_plus3pct_null','lsf_plus3pct_alternative'),('omit_blended','omit_blended_null','omit_blended_alternative'),('opacity_zero','opacity_zero_null','opacity_zero_alternative'),('omit_saturated','omit_saturated_null','omit_saturated_alternative'),('wider_bounds','wider_bounds_null','wider_bounds_alternative')]

def select(name):
 names=[name,name+'_polished']
 if name=='null':names+=['null_cross']
 if name=='alternative':names+=['alternative_from_published']
 results=[json.loads((OUT/f'{n}.json').read_text()) for n in names if (OUT/f'{n}.json').exists()]
 if not results:return None
 converged=[r for r in results if r['optimizer_success']]
 return min(converged or results,key=lambda r:r['chi2'])

def main():
 rows=[]
 for case,an,bn in PAIRS:
  if not (OUT/f'{an}.json').exists() or not (OUT/f'{bn}.json').exists():continue
  a=select(an);b=select(bn);d=a['chi2']-b['chi2'];k=b['npar']-a['npar']
  rows.append(dict(case=case,null_file=a['name']+'.json',alternative_file=b['name']+'.json',n_pixels=a['ndata'],null_parameters=a['npar'],added_parameters=k,null_chi2=a['chi2'],alternative_chi2=b['chi2'],delta_chi2=d,null_chi2_per_nominal_ndf=a['chi2_per_nominal_ndf'],alternative_chi2_per_nominal_ndf=b['chi2_per_nominal_ndf'],nominal_wilks_p=float(chi2.sf(max(0,d),k)),delta_aic_heuristic=2*k-d,delta_bic_heuristic=k*np.log(a['ndata'])-d,null_success=a['optimizer_success'],alternative_success=b['optimizer_success'],null_active_bounds=len(a['active_bounds']),alternative_active_bounds=len(b['active_bounds']),null_optimality=a['optimality'],alternative_optimality=b['optimality']))
 if rows:
  with (OUT/'comparison.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 result=dict(comparisons=rows,interpretation='Paired exploratory one-absorber FeII fits, not a population or time-law test. Nominal Wilks p and AIC/BIC are heuristic: local optimization, active nuisance bounds, component selection and unknown covariance limit calibration.')
 (OUT/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
 if not (OUT/'null.npz').exists():return
 ja=select('null');a=np.load(OUT/(ja['name']+'.npz'))
 plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
 fig,axs=plt.subplots(3,2,figsize=(11,8.8),constrained_layout=True)
 masksummary=[]
 for ax,k in zip(axs.flat,[2260,2344,2374,2382,2586,2600]):
  v=a[f'{k}_v'];f=a[f'{k}_flux'];m=a[f'{k}_model'];e=a[f'{k}_error'];good=a[f'{k}_good'];r=(f[good]-m[good])/e[good]
  ax.step(v,np.where(good,f,np.nan),where='mid',lw=.65,color='0.25',label='Observed ESPRESSO flux')
  ax.plot(v,np.where(good,m,np.nan),color='#bc3f36',lw=1.1,label='Conventional shared-gas fit')
  ax.scatter(v[~good],np.clip(f[~good],-.03,1.08),s=5,c='0.75',label='Author-masked pixels',zorder=1)
  ax.axhline(-.25,color='0.6',lw=.6);ax.fill_between(v,-.30,-.20,color='#78a98a',alpha=.20)
  rr=np.full(len(v),np.nan);rr[good]=-.25+.05*r
  ax.plot(v,rr,lw=.55,color='#285f8f',label='Residuals: 0.05 per sigma')
  per=next(t for t in ja['per_line'] if t['line']==k)
  ax.text(.03,.93,f'Fe II {k}  |  N={good.sum()}  |  $\\chi^2/N$={per["chi2_per_pixel"]:.2f}',transform=ax.transAxes)
  ax.set(ylim=(-.47,1.18),xlabel='Velocity relative to z=1.150793 (km/s)',ylabel='Normalized flux + residual strip')
  adjacent=np.diff(np.flatnonzero(good))==1
  corr=float(np.corrcoef(r[:-1][adjacent],r[1:][adjacent])[0,1]) if adjacent.sum()>2 else None
  masksummary.append(dict(line=k,ndata=int(good.sum()),normalized_residual_mean=float(r.mean()),rms=float(np.sqrt(np.mean(r*r))),adjacent_residual_correlation=corr,interpretation='Residual correlation diagnostic, not an independently measured noise covariance'))
 axs.flat[0].legend(loc='lower left',fontsize=7,framealpha=.95)
 fig.suptitle('Conventional absorption physics on real ESPRESSO data\nFixed laboratory relations; gas and continuum re-fitted; no cosmic-time parameter',fontsize=13)
 for ext in ['png','pdf']:fig.savefig(FIG/f'espresso_conventional_null.{ext}',dpi=180)
 plt.close(fig)
 (OUT/'residual_diagnostics.json').write_text(json.dumps(masksummary,indent=2)+'\n')
 if rows:
  fig,ax=plt.subplots(figsize=(9,4.4),constrained_layout=True)
  names=[r['case'].replace('_',' ') for r in rows];xx=np.arange(len(rows))
  ax.barh(xx,[r['delta_chi2'] for r in rows],color='#386e88',alpha=.9,label='Observed improvement')
  for i,r in enumerate(rows):ax.text(r['delta_chi2']+.7,i,f'+{r["added_parameters"]} shifts',va='center',fontsize=8)
  ax.set(yticks=xx,yticklabels=names,xlabel='Chi-square reduction after adding relative line shifts',title='One-absorber diagnostic: conventional null vs free line shifts')
  ax.set_xlim(0,max(r['delta_chi2'] for r in rows)*1.22)
  ax.legend(loc='best',fontsize=8);ax.invert_yaxis()
  for ext in ['png','pdf']:fig.savefig(FIG/f'espresso_null_robustness.{ext}',dpi=180)
  plt.close(fig)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
