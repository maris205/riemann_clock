#!/usr/bin/env python3
"""Plot empirical controls, clearly separating data from assumed covariance."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'results/noise_covariance'
c=json.loads((OUT/'controls.json').read_text())
s=json.loads((OUT/'local_shift_covariance_sensitivity.json').read_text())
p=c['primary'];acf=np.array(p['acf']);ci=np.array(p['interval_cluster_bootstrap']['acf_95_percentile_intervals'])
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(13.4,4.3),gridspec_kw={'width_ratios':[1.1,1.,1.35]})
ax=axes[0];lags=np.arange(1,11)
ax.errorbar(lags,acf[1:],yerr=[acf[1:]-ci[1:,0],ci[1:,1]-acf[1:]],fmt='o',ms=4,color='#17668b',capsize=2,label='Control intervals (cluster bootstrap)')
ax.plot(lags,acf[1]**lags,'--',color='#b86627',label='AR(1), matched lag 1')
ax.axhline(0,color='.5',lw=.7);ax.set(xlabel='Pixel lag (0.4 km/s per pixel)',ylabel='Autocorrelation',title='Observed nearby controls',xlim=(.5,10.5));ax.legend(fontsize=8)
ax=axes[1]
keys=list(c['by_transition']);x=np.arange(len(keys));rho=[c['by_transition'][key]['acf'][1] for key in keys]
ax.scatter(x,rho,s=35,color='#17668b');ax.axhline(acf[1],ls='--',color='.5',lw=1)
for i,key in enumerate(keys):ax.annotate(f'n={c["by_transition"][key]["interval_count"]}',(i,rho[i]),xytext=(0,7),textcoords='offset points',ha='center',fontsize=8)
ax.set(xticks=x,xticklabels=keys,xlabel='Nearby Fe II transition (Å)',ylabel='Lag-1 correlation',title='Disjoint wavelength controls',ylim=(.2,.55))
ax=axes[2]
names=['diagonal_nominal','ar1_nominal','ma1_nominal','tapered_acf10_nominal']
labels=['Diagonal','AR(1)','MA(1)','Tapered ACF']
rows={r['name']:r for r in s['cases'] if r['relative_singular_value_cutoff']==1e-10}
y=np.arange(4)
ax.scatter([rows[k]['local_score_delta_chi2'] for k in names],y,label='File error scale',s=42,color='#17668b')
scaled=[np.nan]+[rows[k.replace('_nominal','_control_scaled')]['local_score_delta_chi2'] for k in names[1:]]
ax.scatter(scaled,y,marker='s',label='Transferred continuum variance',s=35,color='#b86627')
ax.set(yticks=y,yticklabels=labels,xlabel='Local five-shift score',title='Fixed-baseline covariance sensitivity',xlim=(0,65),ylim=(4.25,-.5));ax.legend(fontsize=8,loc='lower right')
fig.suptitle('25 control intervals / 13,750 pixels: correlation is present outside absorption',fontsize=13)
fig.text(.5,.01,'Exploratory control selection. Right panel is a local nuisance projection, not a nonlinear refit or a discovery significance.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.045,1,.93])
for extension in ['png','pdf']:fig.savefig(OUT/f'controls_and_covariance.{extension}',dpi=180,bbox_inches='tight')
