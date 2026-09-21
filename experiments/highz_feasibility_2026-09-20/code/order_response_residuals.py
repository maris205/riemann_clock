#!/usr/bin/env python3
"""Frozen exploratory diagnostics of native-exposure adjacent-order offsets.
No new line-spread-function model; all uncertainty statements are conditional.
"""
from pathlib import Path
import argparse, json, hashlib, datetime, sys
import numpy as np
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import exposure_analysis as a
ROOT=a.ROOT; OUT=ROOT/'results/order_response'; OUT.mkdir(parents=True,exist_ok=True)
REPORT=ROOT/'reports/order_response_residuals_cn.md'
PFILE=OUT/'residual_protocol.json'
LINES=[2382,2600]
SHA=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x): Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def protocol():
 if PFILE.exists(): raise SystemExit('Protocol already frozen; refusing overwrite.')
 inputs=[ROOT/'results/exposures/fit_rows.json',a.PRO/'metadata.json',a.TEMPLATE_PATH,ROOT/'code/exposure_analysis.py']
 meta=json.loads((a.PRO/'metadata.json').read_text())
 inputs += [ROOT/r['array_file'] for r in meta['exposures']]
 obj=dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='Exploratory follow-up protocol frozen after knowledge of the approximately -60 m/s 2600 order contrast, before running these diagnostic calculations; not a blind preregistration.',lines=LINES,order_contrast='lower order index minus upper order index for the SAME transition; no 2374 reference enters the contrast or its variance',source_hashes={str(p.relative_to(ROOT)):SHA(p) for p in inputs},models=['archived lower_order and upper_order, fixed Gaussian LSF','archived lower_order_free_lsf and upper_order_free_lsf as descriptive sensitivity','new per-trace lower-minus-upper, same fixed shared gas template and nominal Gaussian LSF'],primary_covariance='Conditional independent native ERRDATA; paired variance is sum of disjoint lower and upper order shift variances, no residual chi-square scaling. Unknown extraction/template covariance is not estimated.',trends=dict(primary='Each predictor separately with intercept; report every slope and conditional standard error without selecting p values.',predictors=['MJD in Julian years from first exposure','pipeline BERV km/s','continuum-equivalent per-pixel SNR proxy separately for lower and upper order','log lower/upper SNR ratio','fixed v=0 reference detector pixel coordinate separately for lower and upper order'],pixel_phase='Fractional detector coordinate at fixed template v=0, separately for each order and trace. Both sine and cosine enter together; integer native pixel index is not itself a subpixel phase. This single reference does not describe all blended component phases.',epoch='2018 minus 2019-2020 fixed contrast; observing epoch is not cosmic epoch.',collinearity='Report standardized predictor correlation matrix and singular values; no multivariable causal attribution because BERV, epoch and detector position are not independently varied.',leave_one_out='Recompute constant paired mean and all one-predictor slopes after omitting each of the 17 exposures; no omission becomes preferred.'),residual_diagnostics=dict(sign='model minus data in ERRDATA sigma units after archived per-order fit',velocity_edges_km_s=[-90,-40,0,40,80,150],strength_groups=[['deep_core',0,.1],['core',.1,.3],['shoulder',.3,.9],['continuum',.9,1.01]],strength_assignment='Frozen nominal Gaussian-convolved null_cross template at zero shift integrated over each native pixel; neither observed flux nor newly fitted residual sets categories.',phase_edges=[0,.25,.5,.75,1],summary='Counts, mean and RMS of normalized residuals by line/order/trace and bin, plus per-exposure values. No independence-based pixel p values or best-bin claims.'),trace_fits=136,interpretation='Diagnostics can localize sensitivity but cannot determine a unique instrument/extraction mechanism, establish physical frequency drift, or test 1/log^2 cosmic time. All 17 exposures see the same absorber; shared coadd-derived template prevents independent confirmation.')
 dump(PFILE,obj)
 print('Frozen',PFILE)

def wls(y,var,X):
 y=np.asarray(y);var=np.asarray(var);X=np.asarray(X);Q,R=np.linalg.qr(X/np.sqrt(var)[:,None],mode='reduced');beta=np.linalg.solve(R,Q.T@(y/np.sqrt(var)));cov=np.linalg.solve(R,np.eye(R.shape[0]));cov=cov@cov.T;r=(y-X@beta)/np.sqrt(var);q=float(r@r)
 return dict(beta=beta.tolist(),covariance=cov.tolist(),sigma=np.sqrt(np.diag(cov)).tolist(),chi2=q,ndf=len(y)-X.shape[1],conditional_goodness_p=float(chi2.sf(q,len(y)-X.shape[1])))

def summary(y,var):
 r=wls(y,var,np.ones((len(y),1)));return dict(mean_m_s=r['beta'][0],sigma_m_s=r['sigma'][0],heterogeneity_chi2=r['chi2'],heterogeneity_ndf=r['ndf'],conditional_heterogeneity_p=r['conditional_goodness_p'])

def moments(r):
 r=np.asarray(r,dtype=float)
 return dict(npix=len(r),sum=float(r.sum()),sum_squares=float(r@r),mean=float(r.mean()) if len(r) else None,rms=float(np.sqrt(np.mean(r*r))) if len(r) else None)

def fit_select(record,arr,k,t,order,trace):
 rec=dict(record,segments=[s for s in record['segments'] if s['line']==k and s['trace_index']==trace])
 return a.fit_line(rec,arr,k,t,selection=order,free_lsf=False)

def run():
 p=json.loads(PFILE.read_text())
 for f,h in p['source_hashes'].items():
  if SHA(ROOT/f)!=h:raise RuntimeError('Input changed since frozen protocol: '+f)
 rows=json.loads((a.OUT/'fit_rows.json').read_text());meta=json.loads((a.PRO/'metadata.json').read_text());lookup={(r['configuration'],r['exposure_index'],r['line']):r for r in rows};t=a.Template()
 points=[];rowdesc=[];residual_store=[];tracefits=[]
 for rec in meta['exposures']:
  ei=rec['index']
  with np.load(ROOT/rec['array_file']) as arr:
   for k in LINES:
    pf={}
    for selection in ['lower_order','upper_order']:
     old=lookup[(selection,ei,k)];snrs=[];pix=[];phases=[]
     for s in a.choose_segments(rec,k,selection):
      pre=s['prefix'];good=arr[pre+'_good'];v=arr[pre+'_v'][good];err=arr[pre+'_error'][good]
      d=next(x for x in old['rows'] if x['row']==s['row']);res=np.asarray(d['residual']);cont=d['continuum_amplitude']+d['continuum_slope']*(v-v.mean())/100
      snrs.extend((cont/err).tolist());x0=float(np.interp(0,arr[pre+'_v'],arr[pre+'_pixel']));phase=x0%1;pix.append(x0);phases.append(phase)
      profile=t.spline(k,a.REGIONS[k][2])(t.coordinates(arr,pre,k)[good])@t.weights/2
      st=np.digitize(profile,[.1,.3,.9],right=True);ve=np.digitize(v,p['residual_diagnostics']['velocity_edges_km_s'])-1
      rd=dict(exposure_index=ei,line=k,selection=selection,row=s['row'],trace=s['trace_index'],pixel_at_v0=x0,phase_at_v0=phase,snr_continuum_proxy_median=float(np.median(cont/err)),all=moments(res),velocity_bins=[moments(res[ve==i]) for i in range(5)],strength_bins=[moments(res[st==i]) for i in range(4)])
      rowdesc.append(rd);residual_store.append(dict(line=k,selection=selection,trace=s['trace_index'],phase_bin=min(int(phase*4),3),res=res,ve=ve,st=st))
     pf[selection]=dict(snr_proxy=float(np.median(snrs)),pixel_at_v0=float(np.mean(pix)),pixel_at_v0_by_trace=pix,phase_at_v0_by_trace=phases)
    lower=lookup[('lower_order',ei,k)];upper=lookup[('upper_order',ei,k)];lf=lookup[('lower_order_free_lsf',ei,k)];uf=lookup[('upper_order_free_lsf',ei,k)]
    points.append(dict(exposure_index=ei,line=k,date_obs=rec['date_obs'],mjd=rec['mjd_obs'],berv_km_s=rec['berv_km_s'],early_2018=int(rec['date_obs'].startswith('2018')),delta_m_s=1000*(lower['shift_km_s']-upper['shift_km_s']),variance_m2_s2=1e6*(lower['shift_sigma_km_s']**2+upper['shift_sigma_km_s']**2),delta_free_lsf_m_s=1000*(lf['shift_km_s']-uf['shift_km_s']),variance_free_lsf_m2_s2=1e6*(lf['shift_sigma_km_s']**2+uf['shift_sigma_km_s']**2),**pf))
    for tr in [0,1]:
     for order in ['lower_order','upper_order']:
      f=fit_select(rec,arr,k,t,order,tr);tracefits.append(dict(trace=tr,**f))
  print('exposure',ei,'trace fits',len(tracefits),flush=True)
  dump(OUT/'residual_trace_fits.json',tracefits)
 dump(OUT/'residual_points.json',points);dump(OUT/'residual_rows.json',rowdesc)
 results={}
 for k in LINES:
  ps=[r for r in points if r['line']==k];y=np.array([r['delta_m_s'] for r in ps]);var=np.array([r['variance_m2_s2'] for r in ps]);n=len(y);early=np.array([r['early_2018'] for r in ps]);mjd=np.array([r['mjd'] for r in ps]);berv=np.array([r['berv_km_s'] for r in ps]);sx=np.array([[r[o]['snr_proxy'] for o in ['lower_order','upper_order']] for r in ps]);px=np.array([[r[o]['pixel_at_v0'] for o in ['lower_order','upper_order']] for r in ps]);pred={'years_since_first':(mjd-mjd[0])/365.25,'berv_km_s':berv,'snr_lower':sx[:,0],'snr_upper':sx[:,1],'log_snr_ratio':np.log(sx[:,0]/sx[:,1]),'pixel_lower':px[:,0],'pixel_upper':px[:,1]}
  fits={};loo={}
  for name,x in pred.items():
   center=float(np.average(x,weights=1/var));X=np.column_stack([np.ones(n),x-center]);rr=wls(y,var,X);rr.update(predictor_center=center,units='m/s per predictor unit',x_values=x.tolist());fits[name]=rr;lr=[]
   for j in range(n):
    keep=np.arange(n)!=j;rrj=wls(y[keep],var[keep],X[keep]);lr.append(dict(omitted_exposure=j,slope=rrj['beta'][1],slope_sigma=rrj['sigma'][1]))
   loo[name]=lr
  phasefit={}
  for oi,o in enumerate(['lower_order','upper_order']):
   for tr in [0,1]:
    ph=np.array([r[o]['phase_at_v0_by_trace'][tr] for r in ps]);X=np.column_stack([np.ones(n),np.sin(2*np.pi*ph),np.cos(2*np.pi*ph)]);rr=wls(y,var,X);rr.update(phases=ph.tolist());phasefit[o+'_trace'+str(tr)]=rr
  epoch=wls(y,var,np.column_stack([np.ones(n),early]));epoch['definition']='beta[0] = 2019-2020; beta[1] = 2018 minus 2019-2020'
  corr_names=['years_since_first','early_2018']+list(pred)[1:];xx=np.column_stack([pred['years_since_first'],early]+[pred[j] for j in list(pred)[1:]]);zz=(xx-xx.mean(axis=0))/xx.std(axis=0);corr=np.corrcoef(xx.T);sing=np.linalg.svd(zz,compute_uv=False)
  loomean=[dict(omitted_exposure=j,**summary(np.delete(y,j),np.delete(var,j))) for j in range(n)]
  tracepoints=[];tracesummary={}
  for tr in [0,1]:
   for ei in range(n):
    lo=next(r for r in tracefits if r['line']==k and r['trace']==tr and r['exposure_index']==ei and r['selection']=='lower_order');up=next(r for r in tracefits if r['line']==k and r['trace']==tr and r['exposure_index']==ei and r['selection']=='upper_order')
    tracepoints.append(dict(exposure_index=ei,trace=tr,delta_m_s=1000*(lo['shift_km_s']-up['shift_km_s']),variance_m2_s2=1e6*(lo['shift_sigma_km_s']**2+up['shift_sigma_km_s']**2)))
   tp=[r for r in tracepoints if r['trace']==tr];tracesummary[str(tr)]=summary([r['delta_m_s'] for r in tp],[r['variance_m2_s2'] for r in tp])
  dif=np.array([tracepoints[i]['delta_m_s']-tracepoints[i+n]['delta_m_s'] for i in range(n)]);dv=np.array([tracepoints[i]['variance_m2_s2']+tracepoints[i+n]['variance_m2_s2'] for i in range(n)])
  results[str(k)]=dict(primary=summary(y,var),free_lsf=summary([r['delta_free_lsf_m_s'] for r in ps],[r['variance_free_lsf_m2_s2'] for r in ps]),trends=fits,epoch=epoch,leave_one_out=loomean,leave_one_out_slopes=loo,pixel_phase_trends=phasefit,predictor_correlation=dict(names=corr_names,matrix=corr.tolist(),standardized_singular_values=sing.tolist(),condition_number=float(sing[0]/sing[-1])),trace_points=tracepoints,trace_summary=tracesummary,trace0_minus_trace1=summary(dif,dv))
 aggregate=[]
 for k in LINES:
  for o in ['lower_order','upper_order']:
   for tr in ['all',0,1]:
    ss=[x for x in residual_store if x['line']==k and x['selection']==o and (tr=='all' or x['trace']==tr)]
    rr=np.concatenate([x['res'] for x in ss]);ve=np.concatenate([x['ve'] for x in ss]);st=np.concatenate([x['st'] for x in ss]);pb=np.concatenate([np.full(len(x['res']),x['phase_bin']) for x in ss]);aggregate.append(dict(line=k,selection=o,trace=tr,all=moments(rr),velocity_bins=[moments(rr[ve==i]) for i in range(5)],strength_bins=[moments(rr[st==i]) for i in range(4)],reference_phase_bins=[moments(rr[pb==i]) for i in range(4)]))
 outputs=dict(protocol_sha256=SHA(PFILE),script_sha256=SHA(__file__),source_hashes=p['source_hashes'],results=results,aggregate_residuals=aggregate,trace_fit_count=len(tracefits),trace_fit_successes=sum(r['optimizer_success'] for r in tracefits),trace_bound_active_count=sum(r['bound_active'] for r in tracefits),limitations=p['interpretation'],pixel_phase_note=p['trends']['pixel_phase'])
 dump(OUT/'residual_summary.json',outputs);plot(points,outputs);report(outputs,points)
 print(json.dumps({k:dict(primary=v['primary'],epoch=v['epoch'],trace=v['trace_summary']) for k,v in results.items()},indent=2))

def plot(points,out):
 fig,ax=plt.subplots(2,2,figsize=(12,8));colors={2382:'#3276a5',2600:'#b55a27'}
 for k in LINES:
  ps=[r for r in points if r['line']==k];r=out['results'][str(k)];ys=np.array([q['delta_m_s'] for q in ps]);es=np.sqrt([q['variance_m2_s2'] for q in ps]);c=colors[k]
  ax[0,0].errorbar(np.arange(17)+(0.1 if k==2600 else -.1),ys,yerr=es,fmt='o',ms=3,capsize=2,color=c,label=str(k))
  ax[0,0].axhline(r['primary']['mean_m_s'],color=c,ls='--',lw=.8)
  x=np.array([q['berv_km_s'] for q in ps]);ax[0,1].errorbar(x,ys,yerr=es,fmt='o',ms=3,capsize=2,color=c,label=str(k));b=r['trends']['berv_km_s'];xx=np.linspace(x.min(),x.max(),50);ax[0,1].plot(xx,b['beta'][0]+b['beta'][1]*(xx-b['predictor_center']),color=c,lw=.8)
  loo=r['leave_one_out'];ax[1,0].errorbar(np.arange(17)+(0.1 if k==2600 else -.1),[q['mean_m_s'] for q in loo],yerr=[q['sigma_m_s'] for q in loo],fmt='o',ms=3,capsize=2,color=c,label=str(k))
 for o,marker in [('lower_order','o'),('upper_order','s')]:
  a0=next(x for x in out['aggregate_residuals'] if x['line']==2600 and x['selection']==o and x['trace']=='all');ax[1,1].plot(np.arange(4),[b['rms'] for b in a0['strength_bins']],marker=marker,label=o.replace('_',' '))
 ax[0,0].set(xlabel='Exposure index',ylabel='Lower − upper order [m/s]',title='Same-transition paired offsets');ax[0,0].axvline(8.5,color='gray',ls=':',lw=.8)
 ax[0,1].set(xlabel='Pipeline BERV [km/s]',ylabel='Lower − upper order [m/s]',title='BERV trend is confounded with epoch / pixel position')
 ax[1,0].set(xlabel='Omitted exposure index',ylabel='Mean paired offset [m/s]',title='Leave-one-exposure-out stability')
 ax[1,1].set(xticks=np.arange(4),xticklabels=['F ≤ 0.1','0.1 < F ≤ 0.3','0.3 < F ≤ 0.9','F > 0.9'],ylabel='RMS residual / ERRDATA',title='2600 residuals by frozen template strength')
 for a0 in ax.flat:a0.grid(alpha=.15);a0.legend(fontsize=8)
 fig.suptitle('Adjacent-order diagnostics in 17 native ESPRESSO exposures',fontsize=14);fig.text(.5,.012,'Fixed coadd-derived template; conditional diagonal errors. Observing dates do not sample cosmic evolution.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.035,1,.96]);fig.savefig(OUT/'residual_diagnostics.pdf');fig.savefig(OUT/'residual_diagnostics.png',dpi=180);plt.close(fig)

def report(out,points):
 lines=['# 同一跃迁相邻级次差异：固定诊断协议与残差结构','', '这是在已知 Fe II 2600 存在约 −60 m/s 级次差之后冻结的探索性诊断，不是盲检或新物理检验。所有误差均条件于固定的合并谱气体模板、独立原生 ERRDATA 和指定 Gaussian 仪器响应；未知提取协方差与模板不确定性未计入。','', '本轮能支持的结论是：2600 的级次差不是由单次曝光或仅一条迹线主导；当前有限样本中，未见日期、BERV、信噪代理或参考探测器位置的简单线性趋势能够解释它。强线核心的标准化残差较大，提示应继续检查线形与误差描述，但并不区分仪器、提取、模板或其他常规机制。','', '协议先保存为 `results/order_response/residual_protocol.json`，再运行计算；源数据、旧拟合、模板与代码哈希均记录。此处直接比较同一跃迁的两个级次，2374 参考线完全抵消，不重复加入其方差。','', '| 线 | 低索引−高索引级次 (m/s) | 2018−2019/20 级次差变化 (m/s) | BERV 斜率 (m/s)/(km/s) | 留一曝光均值范围 (m/s) |','|---|---:|---:|---:|---:|']
 for k in LINES:
  r=out['results'][str(k)];s=r['primary'];e=r['epoch'];b=r['trends']['berv_km_s'];lo=[x['mean_m_s'] for x in r['leave_one_out']];lines.append(f"| {k} | {s['mean_m_s']:.3f} ± {s['sigma_m_s']:.3f} | {e['beta'][1]:.3f} ± {e['sigma'][1]:.3f} | {b['beta'][1]:.3f} ± {b['sigma'][1]:.3f} | {min(lo):.3f} 至 {max(lo):.3f} |")
 lines += ['', '## 两条迹线独立拟合同一跃迁的相邻级次','',f"共 {out['trace_fit_count']} 次新增固定模板、固定 Gaussian 宽度拟合，{out['trace_fit_successes']} 次正常终止，{out['trace_bound_active_count']} 次碰速度边界。碰界条目是短曝光 index=10 的 2382、低索引级次、迹线 0，速度 +1.4999999 km/s；该点保留，局部误差不是无约束 Gaussian 置信区间，2600 的迹线拟合没有碰界。两迹线原生样本不重叠；共享模板及系统误差仍使其不是独立天体物理确认。",'', '| 线 | 迹线 0 (m/s) | 迹线 1 (m/s) | 配对迹线 0−1 (m/s) |','|---|---:|---:|---:|']
 for k in LINES:
  r=out['results'][str(k)];ss=[r['trace_summary'][str(t)] for t in [0,1]]+[r['trace0_minus_trace1']];lines.append('| '+str(k)+' | '+' | '.join(f"{s['mean_m_s']:.3f} ± {s['sigma_m_s']:.3f}" for s in ss)+' |')
 lines += ['', '配对迹线差按每次曝光的迹线差及其方差合并，权重不同于两条迹线各自合并的权重，因此表中最后一列不必恰好等于前两列均值相减。','', '## 预先列出的全部单预测量趋势','', '每一项单独加截距拟合，并保留所有结果；不挑最小 p 值，也不将单变量相关解释为因果。SNR 是拟合连续谱/ERRDATA 的逐像素中位数代理，像素位置是固定 v=0 参考位置在两条迹线上的平均插值坐标。坐标变化主要反映同一谱线在探测器上的移动。','', '| 预测量 | 2382 斜率 ± 条件误差 | 2600 斜率 ± 条件误差 |','|---|---:|---:|']
 for name in out['results']['2600']['trends']:
  rr=[out['results'][str(k)]['trends'][name] for k in LINES];lines.append('| '+name+' | '+' | '.join(f"{r['beta'][1]:.5g} ± {r['sigma'][1]:.5g}" for r in rr)+' |')
 lines += ['', '各斜率单位均为 m/s 每预测量单位；时间单位为 Julian 年。完整协方差、各次留一斜率、每个级次和迹线的正弦/余弦相位回归均保存在 JSON 中。没有用曝光年份代替吸收系统的宇宙年龄。','', '## 混杂、像素相位和残差结构','']
 for k in LINES:
  r=out['results'][str(k)];cc=r['predictor_correlation'];names=cc['names'];c=np.asarray(cc['matrix']);lines.append(f"- {k}: BERV 与低/高级次像素位置的 Pearson 相关分别为 {c[names.index('berv_km_s'),names.index('pixel_lower')]:.7f}、{c[names.index('berv_km_s'),names.index('pixel_upper')]:.7f}；BERV 与年份分组相关 {c[names.index('berv_km_s'),names.index('early_2018')]:.5f}；标准化设计条件数 {cc['condition_number']:.3g}。")
 lines += ['', '像素相位严格定义为固定 v=0 参考点的插值探测器坐标的小数部分，而不是原生整数像素编号的余数。每个级次、迹线分别记录；相位是环变量，因此正弦和余弦成对拟合。相位回归的因变量仍是两条迹线合并拟合得到的相邻级次差，分别使用各迹线的参考相位作为预测量；不是迹线专属位移的相位回归。单个参考点不能代表全部混合气体分量的相位。相位分箱残差仅展示样本组成敏感性，不能诊断独立的像素内响应机制。','', '归一化残差为 model−data；按预设速度段以及固定模型通量阈值 0.1、0.3、0.9 分组，不按观测噪声或残差选区。这里不对大量像素采用独立性假设构造发现 p 值。','', '| 线/级次 | F≤0.1 RMS | 0.1<F≤0.3 RMS | 0.3<F≤0.9 RMS | F>0.9 RMS |','|---|---:|---:|---:|---:|']
 for rr in out['aggregate_residuals']:
  if rr['trace']=='all':lines.append('| '+str(rr['line'])+'/'+rr['selection']+' | '+' | '.join(f"{b['rms']:.4f} (n={b['npix']})" if b['npix'] else '无样本' for b in rr['strength_bins'])+' |')
 lines += ['', '完整速度分箱、强度分箱、参考相位分箱及逐曝光逐行统计保存在 `residual_summary.json` 和 `residual_rows.json`。这些残差描述不检验未知相关噪声的频率结构，也不能单凭 RMS 大小在气体模板误差、标定、提取、LSF 或像素响应之间作选择。','', '本轮不修改旧主结果、数据掩膜或响应模型，不选择曝光剔除方案；新增迹线拟合仍使用已有原生像素与冻结气体模板。两个级次的波长覆盖、像素采样与保留掩膜也不同，因此级次差本身不能单独识别波长标定偏差。跨曝光/迹线一致并不排除所有曝光共有的处理或模板偏差。图为 `results/order_response/residual_diagnostics.pdf` / `.png`。','']
 REPORT.write_text('\n'.join(lines))
def union_intervals(left,right,tolerance=1e-10):
    order=np.argsort(left); merged=[]
    for j in order:
        lo=float(left[j]);hi=float(right[j])
        if merged and lo<=merged[-1][1]+tolerance: merged[-1][1]=max(hi,merged[-1][1])
        else: merged.append([lo,hi])
    return np.asarray(merged,dtype=float).reshape(-1,2)

def freeze_common():
    path=OUT/'residual_common_support_protocol.json'
    if path.exists(): raise SystemExit('Common-support protocol already frozen.')
    prior=json.loads(PFILE.read_text())
    obj=dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='Exploratory sensitivity requested after original order-trend diagnostics; frozen before common-support masks/fits were computed.',source_hashes=prior['source_hashes'],lines=LINES,number_of_fits=68,mask_definition='For each exposure and same transition, union all originally valid native wavelength pixel-bin intervals in each order across both traces. Retain an originally good pixel only if its whole wavelength bin lies inside a merged valid-support interval of the opposite order. Single pass using original valid support; no iterative erosion across mismatched native grids. Every retained bin lies in the intersection of original valid supports, although final native-bin unions can have small differing edge footprints.',interval_boundary_tolerance_A=1e-10,fit='Same frozen null_cross template, nominal Gaussian widths, native ERRDATA and separately profiled amplitude/slope/zero for each row as original per-order diagnostic.',scope='Two orders x two transitions x17 exposures. No 2374 reference. All valid flux/errors/masks and original outputs retained unchanged; temporary arrays only.',comparison='Each common-support same-line order difference has variance equal to sum of its disjoint order variances. Change relative to original estimator is descriptive, without independent-estimator error combination because photons overlap.',selection='No observed flux/residual/significance-based masking or exposure exclusion; whole-bin geometry only.',interpretation='Tests differing support/mask sensitivity, does not equalize pixel sampling, per-pixel SNR, instrumental response or noise covariance; not a calibration-only test.')
    dump(path,obj);print('Frozen',path)

def common_support():
    path=OUT/'residual_common_support_protocol.json';protocol=json.loads(path.read_text())
    for f,h in protocol['source_hashes'].items():
        if SHA(ROOT/f)!=h: raise RuntimeError('Common-support source changed: '+f)
    meta=json.loads((a.PRO/'metadata.json').read_text());template=a.Template();fits=[];masks=[];points=[]
    for record in meta['exposures']:
        with np.load(ROOT/record['array_file']) as raw: arrays={k:raw[k].copy() for k in raw.files}
        for k in LINES:
            orders={order:a.choose_segments(record,k,order) for order in ['lower_order','upper_order']};support={}
            for order,segments in orders.items():
                left=np.concatenate([arrays[s['prefix']+'_left'][arrays[s['prefix']+'_good']] for s in segments]);right=np.concatenate([arrays[s['prefix']+'_right'][arrays[s['prefix']+'_good']] for s in segments]);support[order]=union_intervals(left,right,protocol['interval_boundary_tolerance_A'])
            for order,segments in orders.items():
                opposite=support['upper_order' if order=='lower_order' else 'lower_order'];tol=protocol['interval_boundary_tolerance_A']
                for seg in segments:
                    pre=seg['prefix'];left=arrays[pre+'_left'];right=arrays[pre+'_right'];before=arrays[pre+'_good'].copy();inside=np.zeros(len(left),bool)
                    for lo,hi in opposite: inside|=(left>=lo-tol)&(right<=hi+tol)
                    arrays[pre+'_good']=before&inside
                    masks.append(dict(exposure_index=record['index'],line=k,selection=order,row=seg['row'],prefix=pre,original_good=int(before.sum()),retained_good=int(arrays[pre+'_good'].sum()),original_indices_kept=np.flatnonzero(arrays[pre+'_good']).tolist(),opposite_original_support_A=opposite.tolist()))
            fs=[]
            for order in ['lower_order','upper_order']:
                fit=a.fit_line(record,arrays,k,template,selection=order,free_lsf=False);fs.append(fit);fits.append(fit)
            points.append(dict(exposure_index=record['index'],date_obs=record['date_obs'],line=k,delta_m_s=1000*(fs[0]['shift_km_s']-fs[1]['shift_km_s']),variance_m2_s2=1e6*(fs[0]['shift_sigma_km_s']**2+fs[1]['shift_sigma_km_s']**2)))
        print('common support exposure',record['index'],flush=True)
    old=json.loads((OUT/'residual_summary.json').read_text());result={}
    for k in LINES:
        pts=[x for x in points if x['line']==k];rr=summary([x['delta_m_s'] for x in pts],[x['variance_m2_s2'] for x in pts]);rr.update(original_mean_m_s=old['results'][str(k)]['primary']['mean_m_s'],descriptive_change_from_original_m_s=rr['mean_m_s']-old['results'][str(k)]['primary']['mean_m_s']);result[str(k)]=rr
    output=dict(protocol_sha256=SHA(path),script_sha256=SHA(__file__),points=points,summary=result,total_fits=len(fits),successes=sum(x['optimizer_success'] for x in fits),bound_active_count=sum(x['bound_active'] for x in fits),total_original_good=sum(x['original_good'] for x in masks),total_retained_good=sum(x['retained_good'] for x in masks),retained_by_line={str(k):dict(original=sum(x['original_good'] for x in masks if x['line']==k),retained=sum(x['retained_good'] for x in masks if x['line']==k)) for k in LINES},notes=protocol['interpretation'],comparison=protocol['comparison'])
    dump(OUT/'residual_common_support_masks.json',masks);dump(OUT/'residual_common_support_fits.json',fits);dump(OUT/'residual_common_support_summary.json',output)
    # Preserve base report generated from its original numerical result and append this control.
    report(old,json.loads((OUT/'residual_points.json').read_text()))
    lines=['','## 附加的共同有效波长支持控制','','本项在查看原级次差及趋势诊断后另行冻结协议（`residual_common_support_protocol.json`），不属于盲检。每个曝光、跃迁先合并某一级次两条迹线的原有有效像素区间，保留另一侧完全落在该有效支持中的原生像素整段；两侧均使用原掩膜的一次性支持，不按通量选择，也不迭代侵蚀不同采样网格。所有旧数组和结果不变。','',f"共 {len(fits)} 次拟合，{output['successes']} 次正常终止、{output['bound_active_count']} 次碰速度边界。两条线原有有效样本 {output['total_original_good']} 个，保留 {output['total_retained_good']} 个。这里只使保留的像素落在共同原始有效支持内；最终原生像素边缘覆盖仍可能略不同，信噪、采样及响应未被等同化。",'', '| 线 | 原级次差 (m/s) | 共同支持级次差 (m/s) | 均值变化，描述性 (m/s) |','|---|---:|---:|---:|']
    for k in LINES:
        r=result[str(k)];lines.append(f"| {k} | {r['original_mean_m_s']:.3f} | {r['mean_m_s']:.3f} ± {r['sigma_m_s']:.3f} | {r['descriptive_change_from_original_m_s']:.3f} |")
    lines += ['', '每行的共同支持误差仅来自两个互不重叠级次的条件方差之和。它与原估计使用重叠光子，所以不为“新−旧”均值变化提供独立误差相加的检验。即使共同支持后位移仍存在，也不能据此排除残余标定、采样、线形、共同模板或误差相关因素。精确保留索引及对侧原有支持区间均保存在 `residual_common_support_masks.json`。','']
    REPORT.write_text(REPORT.read_text()+'\n'.join(lines));old['script_sha256']=SHA(__file__);dump(OUT/'residual_summary.json',old)
    print(json.dumps(output,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--freeze-protocol',action='store_true');ap.add_argument('--freeze-common-support',action='store_true');ap.add_argument('--common-support',action='store_true');args=ap.parse_args()
    if args.freeze_protocol: protocol()
    elif args.freeze_common_support: freeze_common()
    elif args.common_support: common_support()
    else: run()
