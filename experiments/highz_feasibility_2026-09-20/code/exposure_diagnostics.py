#!/usr/bin/env python3
"""Controls and displays for the archived native-exposure fixed-template check."""
from pathlib import Path
import json,sys
import numpy as np
from scipy.stats import chi2
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import exposure_analysis as a
ROOT=a.ROOT;OUT=a.OUT

def paired_diagnostic(rows,first,second,anchor_cancels=False):
 results=[]
 for ei in range(17):
  q={cfg:{r['line']:r for r in rows if r['configuration']==cfg and r['exposure_index']==ei} for cfg in [first,second]}
  if anchor_cancels:
   y=np.array([(q[first][k]['shift_km_s']-q[second][k]['shift_km_s'])*1000 for k in [2382,2600]])
   cov=np.diag([(q[first][k]['shift_sigma_km_s']**2+q[second][k]['shift_sigma_km_s']**2)*1e6 for k in [2382,2600]])
  else:
   p1=a.relative(list(q[first].values()));p2=a.relative(list(q[second].values()));y=np.array(p1['relative_m_s'])-p2['relative_m_s'];cov=np.array(p1['covariance_m2_s2'])+p2['covariance_m2_s2']
  results.append(dict(exposure_index=ei,relative_m_s=y.tolist(),covariance_m2_s2=cov.tolist()))
 out=a.combine(results);mean=np.array(out['mean_relative_m_s']);cov=np.array(out['conditional_covariance_m2_s2']);q=float(mean@np.linalg.solve(cov,mean));out.update(points=results,zero_difference_chi2=q,zero_difference_ndf=2,conditional_zero_difference_p=float(chi2.sf(q,2)),definition=first+' minus '+second,shared_anchor_cancelled=anchor_cancels)
 return out

def main():
 rows=json.loads((OUT/'fit_rows.json').read_text());summary=json.loads((OUT/'summary.json').read_text());meta=json.loads((a.PRO/'metadata.json').read_text());template=a.Template()
 # Apply exactly the same fixed-shape/continuum method to the previously used coadd.
 arr={};segments=[]
 for k in a.KEYS:
  ref=np.average(a.ATOMS[k][:,0],weights=a.ATOMS[k][:,1]);lo,hi,_=a.REGIONS[k];ix=np.flatnonzero((a.COADD_WAVE>=lo)&(a.COADD_WAVE<=hi));pre=f'{k}_coadd';wave=a.COADD_WAVE[ix]
  arr[pre+'_wave']=wave;arr[pre+'_left']=(wave+a.COADD_WAVE[ix-1])/2;arr[pre+'_right']=(wave+a.COADD_WAVE[ix+1])/2;arr[pre+'_v']=a.C*np.log(wave/(ref*(1+a.Z)));arr[pre+'_flux']=a.COADD_DATA[0,ix];arr[pre+'_error']=a.COADD_DATA[1,ix];arr[pre+'_good']=(a.COADD_DATA[4,ix]==1)&(a.COADD_DATA[1,ix]>0)
  segments.append(dict(prefix=pre,line=k,row=0,order_index=0,trace_index=0))
 record=dict(index=-1,date_obs='combined',segments=segments);coaddrows=[a.fit_line(record,arr,k,template) for k in a.KEYS];coaddrelative=a.relative(coaddrows)
 paired_trace=paired_diagnostic(rows,'trace0','trace1');paired_order=paired_diagnostic(rows,'lower_order','upper_order',True);paired_order_lsf=paired_diagnostic(rows,'lower_order_free_lsf','upper_order_free_lsf',True)
 bounds=[dict(configuration=r['configuration'],exposure_index=r['exposure_index'],line=r['line'],fwhm_km_s=r['fwhm_km_s']) for r in rows if r['bound_active']]
 out=dict(coadd_same_fixed_template=coaddrelative,coadd_same_fixed_template_rows=coaddrows,trace_difference=paired_trace,order_difference=paired_order,order_difference_free_lsf=paired_order_lsf,bound_active_fits=bounds,total_good_native_pixels=sum(s['ngood'] for e in meta['exposures'] for s in e['segments']),total_published_mask_transferred_rejections=sum(s['additional_upl_masked_pixels'] for e in meta['exposures'] for s in e['segments']),total_atmospheric_mask_transferred_rejections=sum(s['additional_atmospheric_masked_pixels'] for e in meta['exposures'] for s in e['segments']),conditional_covariance_warning='All reported error bars and p-values condition on the fixed coadd-derived gas template and diagonal native ERRDATA; order and trace comparisons are diagnostic, post-inspection, and not calibrated discovery tests.')
 (OUT/'diagnostics.json').write_text(json.dumps(out,indent=2)+'\n')
 fig,axs=plt.subplots(2,2,figsize=(12.4,8.1),gridspec_kw={'width_ratios':[1.45,1]})
 points=summary['primary']['points'];x=np.arange(17)
 for j,k in enumerate([2382,2600]):
  ax=axs[j,0];y=np.array([p['relative_m_s'][j] for p in points]);err=np.sqrt([p['covariance_m2_s2'][j][j] for p in points]);ax.errorbar(x,y,yerr=err,fmt='o',ms=4,capsize=2,color=['#286c9b','#9b4f27'][j]);ax.axhline(0,color='gray',lw=.7);ax.axvline(8.5,color='gray',ls=':');ax.axhline(summary['primary']['all_exposures']['mean_relative_m_s'][j],color='black',lw=1,ls='--');ax.set_ylabel(f'Fe II {k} − 2374 [m/s]');ax.set_xticks(x);ax.set_xticklabels([p['date_obs'][5:10] for p in points],rotation=60,fontsize=8);ax.set_title('17 native-pixel exposures' if j==0 else '2018 (left)  |  2019–2020 (right)',fontsize=11);ax.grid(alpha=.15)
  ax=axs[j,1];configs=['primary','free_lsf','trace0','trace1','lower_order','upper_order','published_template'];labels=['All rows, fixed LSF','All rows, free LSF','Trace 0','Trace 1','Lower order','Upper order','Published gas template']
  yy=np.arange(len(configs));means=[summary[c]['all_exposures']['mean_relative_m_s'][j] for c in configs];errs=[summary[c]['all_exposures']['conditional_sigma_m_s'][j] for c in configs];ax.errorbar(means,yy,xerr=errs,fmt='o',capsize=3,color=['#286c9b','#9b4f27'][j]);ax.set_yticks(yy);ax.set_yticklabels(labels,fontsize=9);ax.invert_yaxis();ax.axvline(0,color='gray',lw=.8);ax.set_xlabel('Conditional mean relative offset [m/s]');ax.grid(axis='x',alpha=.15)
 fig.suptitle('ESPRESSO exposure consistency — fixed shared gas templates',fontsize=15);fig.text(.5,.014,'Conditional diagonal-noise errors; template learned from the same coadd. These are consistency checks, not independent detections.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.04,1,.95]);fig.savefig(OUT/'exposure_consistency.png',dpi=180);fig.savefig(OUT/'exposure_consistency.pdf');plt.close(fig)
 print(json.dumps({k:v for k,v in out.items() if k not in ['coadd_same_fixed_template_rows']},indent=2))
if __name__=='__main__':main()
