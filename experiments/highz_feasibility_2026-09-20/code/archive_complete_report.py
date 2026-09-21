#!/usr/bin/env python3
"""Describe all four remaining candidates, including inadequate outcomes."""
import json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from archive_complete_model import ROOT,OUT,Model,read_json,sha

TARGETS=['J053007-250329','J064326-504112','J225719-100104','J233156-090802']
Z={'J053007-250329':2.141,'J064326-504112':2.659,'J225719-100104':1.836,'J233156-090802':2.143}

def target_record(t):
 folder=OUT/t;ordinary=read_json(folder/'ordinary_model_summary.json')
 if t=='J225719-100104' and (folder/'selected_pair.json').exists():
  p=read_json(folder/'selected_pair.json');h0=read_json(folder/(p['null']+'.json'))
 elif (folder/'selected_pair.json').exists():
  p=read_json(folder/'selected_pair.json');h0=read_json(folder/(p['null']+'.json'))
 else:p=None;h0=read_json(folder/(ordinary['selected_null']+'.json'))
 record={'target':t,'z_abs':Z[t],'ndata':h0['ndata'],'ncomp':h0['ncomp'],
   'lines':h0['configuration']['lines'],'window':h0['configuration']['window'],
   'primary_or_recovered_null':h0['name'],'null_chi2':h0['chi2'],'null_nominal_ndf':h0['nominal_ndf'],
   'null_stop_success':h0['optimizer_success'],'null_stationarity':h0['stationarity'],
   'null_active_bounds':h0['active_bounds'],'selected_at_count_search_upper_edge':True,
   'source_json_sha256':sha(folder/(h0['name']+'.json'))}
 if p:
  profile=read_json(folder/'pair_profile.json')
  record.update(disposition='conditional_relative_profile_only',alternative=p['alternative'],
    delta_chi2=p['delta_chi2'],extra_shifts=p['extra_shifts'],D_m_s=p['local_pair_uncertainty']['shift_m_s'],
    conditional_sigma_m_s=p['local_pair_uncertainty']['rank_sensitivity'][1]['conditional_sigma_m_s'],
    profile_grid=profile['grid'],pair_file=str((folder/'selected_pair.json').relative_to(ROOT)),
    pair_sha256=sha(folder/'selected_pair.json'))
 else:record.update(disposition='ordinary_model_not_adequate_for_precision_interpretation',D_m_s=None,
      conditional_sigma_m_s=None,delta_chi2=None,extra_shifts=None)
 # Failure-case H1 is a diagnostic of model inadequacy, not a selected precision datum.
 failure=folder/'failure_pair'
 if failure.exists():
  summaries=[p for p in failure.glob('*.json') if 'summary' in p.name]
  record['failure_case_diagnostic_files']=[str(p.relative_to(ROOT)) for p in summaries]
  record['failure_case_diagnostic_interpretation']='Use separately; no promotion into precision sample and no population-null inference from the adequacy gate.'
 return record

def j2257_figures(record):
 folder=OUT/record['target'];pair=read_json(folder/'selected_pair.json')
 null=read_json(folder/(pair['null']+'.json'));alt=read_json(folder/(pair['alternative']+'.json'))
 h0=np.load(folder/(pair['null']+'.npz'));h1=np.load(folder/(pair['alternative']+'.npz'))
 keys=null['configuration']['lines'];fig,axes=plt.subplots(2,len(keys),figsize=(13,6),sharex='col',gridspec_kw={'height_ratios':[3,1]})
 for j,k in enumerate(keys):
  v=h0[f'{k}_v'];g=h0[f'{k}_good'];f=h0[f'{k}_flux'];e=h0[f'{k}_error']
  axes[0,j].plot(v[g],f[g],color='.3',lw=.65,label='Native data');axes[0,j].fill_between(v[g],f[g]-e[g],f[g]+e[g],color='.7',alpha=.25)
  for data,color,label in [(h0,'#247a80','Shared gas H0'),(h1,'#d88227','Relative shifts H1')]:
   model=data[f'{k}_model'];axes[0,j].plot(v,model,color=color,lw=1.0,label=label);axes[1,j].plot(v[g],(model[g]-f[g])/e[g],color=color,lw=.6)
  axes[0,j].set_title(f'Fe II {k}');axes[0,j].set_ylim(-.08,1.18);axes[1,j].axhline(0,color='.5',lw=.5);axes[1,j].set_ylim(-6,6);axes[1,j].set_xlabel('Velocity (km/s)')
 axes[0,0].legend(fontsize=8);axes[0,0].set_ylabel('Normalized flux');axes[1,0].set_ylabel('Residual / expected error')
 fig.suptitle(f'J225719-100104 z=1.836: 30 shared gas components\nFinite multistart result after retaining a lower capped endpoint; no calibrated frequency or time measurement',fontsize=11)
 fig.tight_layout(rect=[0,0,1,.9]);fig.savefig(folder/'physical_profiles.pdf');fig.savefig(folder/'physical_profiles.png',dpi=150);plt.close(fig)
 profile=read_json(folder/'pair_profile.json');x=[q['value_m_s'] for q in profile['grid']];y=[q['delta_from_free'] for q in profile['grid']]
 fig,ax=plt.subplots(figsize=(6,4));ax.plot(x,y,'-',color='#247a80',lw=.8)
 for stopped,marker,color,label in [(True,'o','#247a80','Objective stopping criterion met'),(False,'x','#b64e36','Evaluation cap reached')]:
  points=[q for q in profile['grid'] if q['optimizer_success']==stopped]
  if points:ax.scatter([q['value_m_s'] for q in points],[q['delta_from_free'] for q in points],marker=marker,color=color,label=label,zorder=3)
 ax.legend(fontsize=8);ax.axhline(1,color='.5',ls='--');ax.axhline(4,color='.5',ls=':')
 ax.set(xlabel='Fe II2382 relative to2374 (m/s)',ylabel='Profile objective minus free endpoint',title='Sparse nuisance-reoptimized diagnostic\nBounds/architecture/calibration limits remain');fig.tight_layout();fig.savefig(folder/'common_pair_profile.pdf');fig.savefig(folder/'common_pair_profile.png',dpi=150);plt.close(fig)
 quad=[]
 for r in [null,alt]:
  scores={}
  for os in [21,55]:
   model=Model(r['configuration'],r['ncomp'],r['free_shifts'],os);res=model.fun(np.asarray(r['parameters']));scores[str(os)]=float(res@res)
  quad.append({'name':r['name'],'chi2_at_fixed_parameters':scores,'delta_55_minus_21':scores['55']-scores['21']})
 (folder/'quadrature_check.json').write_text(json.dumps({'comparisons':quad,'interpretation':'Fixed-parameter integration check, not reoptimization.'},indent=2)+'\n')

def main():
 rows=[target_record(t) for t in TARGETS]
 j2257_figures(next(r for r in rows if r['target']=='J225719-100104'))
 fig,axes=plt.subplots(2,2,figsize=(11,8))
 for ax,row in zip(axes.flat,rows):
  folder=OUT/row['target'];summary=read_json(folder/'ordinary_model_summary.json');points=[]
  for r in summary['by_count']:
   points.append((r['ncomp'],r['chi2']/r['nominal_ndf']))
  if row['target']=='J225719-100104':
   for n in [26,30]:
    data=[read_json(p) for p in folder.glob(f'n{n}_null_extension_s*_os9.json')]
    best=min(data,key=lambda r:r['chi2']);points.append((n,best['chi2_per_ndf']))
  points=sorted(points);ax.plot(*np.asarray(points).T,'o-',color='#247a80')
  ax.axhline(1.5,color='#ce8f49',ls='--',label='Total diagnostic gate');ax.axhline(1,color='.5',ls=':')
  ax.set(yscale='log',xlabel='Shared gas components',ylabel='H0 chi-square / nominal dof',title=f'{row["target"]}\nz_abs={row["z_abs"]:.3f}')
  ax.grid(alpha=.15)
 axes[0,0].legend(fontsize=8);fig.suptitle('All remaining candidates retained, including model inadequacy\nA low residual statistic alone does not establish physical validity',fontsize=12);fig.tight_layout(rect=[0,0,1,.91])
 fig.savefig(OUT/'all_four_gas_complexity.pdf');fig.savefig(OUT/'all_four_gas_complexity.png',dpi=150);plt.close(fig)
 summary={'created_utc':datetime.now(timezone.utc).isoformat(),
   'scope':'All four previously unmodeled automatic candidates; adequacy-gated precision exploration plus separate failed-model H1 diagnostics. This selection cannot establish a population null.',
   'observed_quantity':'D=v(FeII2382)-v(FeII2374), relative-region registration conditional on nuisance/model assumptions; not alpha.',
   'records':rows,
   'previous_two_latest_source':'results/common_pair/harmonized_measurements.json',
   'previous_two_latest_source_sha256':sha(ROOT/'results/common_pair/harmonized_measurements.json'),
   'all_six_disposition_policy':'Retain all six candidates. Two failed-primary ordinary-model candidates are not silently dropped; no physical time-law constraint is claimed.'}
 (OUT/'campaign_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 lines=['# 剩余四个档案候选：有界共同气体模型与失败样本诊断','',
 '**本轮把此前自动筛选出的四个剩余目标全部推进到实际物理谱形建模，保留不充分结果。数据与模型不允许把“未通过常规模型质量门槛”解释成“总体没有频率变化”，也没有构成 α 或对数宇宙时间律的测量。**','',
 '## 范围与固定方法','',
 '全部使用已归档 UVES SQUAD 的原生合并谱、源像素有效掩码以及作者建议用于χ²拟合的 expected-fluctuation 数组。没有按拟合残差剪除像素。窗口及谱线先依据实际谱形、共同吸收结构和污染风险确定；这些仍是探索性选择，不是盲法确认。气体的N、v、b由多条Fe II共同约束；逐线连续谱、零点与Gaussian仪器宽度一起优化。','',
 '新实现调用此前经审查的同位素Voigt光深计算，先合计光深再计算exp(−τ)，随后仪器卷积和原生像素积分。比较 H0/H1 时保留相同像素和误差。连续谱常数界限±0.10、斜率±0.03／100km/s、零点±0.02、仪器FWHM为名义值0.6–1.4倍，气体b为0.5–30km/s。多个分量数和三个起点在固定预算内计算；成功停止与严格平稳、全局最优分开记录。','',
 '用于允许精度诊断的事前质量门槛为总体χ²/名义自由度≤1.5、每条线χ²/像素≤1.8且优化器正常停止。**这个门槛可能排除真实的相对频移，因此不能在通过门槛的子样本上宣称总体不变。** 为避免这个误读，对两个门槛失败目标另授权固定主像素/主模型的有限H1及H0交叉拟合，作为失配解释诊断；它们不会因为Δχ²大就被提升为精密测量。','',
 '## 四个目标的常规基准与有限对比','',
 '|目标|吸收红移|像素|气体分量|H0χ²/ν|主分析状态|','|---|---:|---:|---:|---:|---|']
 for r in rows:lines.append(f'|{r["target"]}|{r["z_abs"]:.3f}|{r["ndata"]}|{r["ncomp"]}|{r["null_chi2"]:.6f}/{r["null_nominal_ndf"]}|'+('仅条件位移诊断' if r['D_m_s'] is not None else '普通基准不足；保留失败样本诊断')+'|')
 lines += ['', '|具备条件诊断的目标|D=v2382−v2374 (m/s)|局部条件σ (m/s)|Δχ² / 新参数|','|---|---:|---:|---:|']
 for r in rows:
  if r['D_m_s'] is not None:lines.append(f'|{r["target"]}|{r["D_m_s"]:+.3f}|{r["conditional_sigma_m_s"]:.3f}|{r["delta_chi2"]:.6f}/{r["extra_shifts"]}|')
 lines += ['',
 '局部σ将2382的Jacobian方向投影到其余全部气体、连续谱、零点、仪器宽度和其他相对位移干扰方向的正交补；同时检验三个SVD截断。它没有包括活动边界、气体架构选择、波长标定、混合线或完整像素协方差，不能当作物理常数的总误差。另保存了真正重新优化干扰参数的稀疏位移剖面。','',
 '## 为什么必须保留普通解释和数值路径','',
 '- **J053007−250329：**加入可见的弱1611以约束饱和，常规基准需要16个分量才通过指定门槛。额外三个位移改善有限；16仍是最初分量网格上限，零点有活动边界，一阶最优性指标也未达到严格平稳。','',
 '- **J064326−504112：**四线主模型的1611残差未过逐线门槛。事前记录了5897Å天空风险，但最大残差位于5893.307Å，不能把失配简单归于那条天空线。另作固定16分量、仅H0的去1611控制后常规残差降低；舍去弱线也损失饱和约束，因此这不是已经证明污染。主四线H1仅作为单独失败案例诊断。','',
 '- **J225719−100104：**原6–22分量网格的AICc持续降低且常规拟合未过门槛，首次位移拟合前限定扩展至26/30分量。在30分量发现一个未完成端点χ²≈422.19，比原先“正常停止”的常规端点χ²≈463.72更低；不能据后者制造额外位移偏好。原端点全部保存，对更低起点进行一次有明确预算的续算后再比较两假设。','',
 '- **J233156−090802：**宽系统的Fe1611会叠入1608谱区红侧，精确间隔约512.033455km/s。模型显式合计同气体的1608和1611光深，并且只保留一份该谱区的似然像素，全部1176个源索引唯一。26个分量后仍未通过总体和2382逐线门槛，因此其H1结果仅是失配诊断。','',
 '## 仍不可作出的推断','',
 '所有四个目标的选择都落在本轮最终分量搜索上端，未证明架构平台。有限多起点会发现不同的局部极小，目标函数停止条件并不等于所有方向已充分平稳。两个较紧凑旧目标的共同谱线对也经过新的有界剖面分析，其中J004131出现更低的常规与替代模型极小；最新数值以 `results/common_pair/harmonized_measurements.json` 和对应报告为准，旧结果作为历史保留。','',
 '六个自动候选均须列在最终样本状态表中。普通模型不足、强线饱和或污染未定的目标，不能既被剔除又被当作支持“没有变化”的证据；反过来，其自由位移收益也不能自动证明新物理。当前数据可用于检验具体谱形模型、评估识别条件与设计后续观测，但尚未提供独立原子响应或每目标差分标定，从而不能将D拟合成物理的1/ln²(t/t*)规律。','',
 '## 报告与复算入口','',
 '- `reports/archive_complete_j0530_results_cn.md`','- `reports/archive_complete_j0643_results.md`','- `reports/archive_complete_j2257_results_cn.md`','- `reports/archive_complete_j2331_results_cn.md`',
 '- `reports/archive_complete_validation.md`：独立原始FITS、像素唯一性、原子不透明度、导数和保存目标复核。',
 '- `results/archive_complete/campaign_summary.json`：所有四个候选与原两个目标最新摘要的来源链接。',
 '- `results/archive_complete/all_four_gas_complexity.pdf/.png`：常规残差随分量数的变化。','',
 '来源：[SQUAD数据论文](https://doi.org/10.1093/mnras/sty2834)、[原子数据](https://doi.org/10.1093/mnras/stt2204)。这些是历史公开观测的新增计算分析，不是2026年新拍摄的光谱。','']
 (ROOT/'reports/archive_complete_results_cn.md').write_text('\n'.join(lines))
 r=next(x for x in rows if x['target']=='J225719-100104')
 targetlines=['# J225719−100104：共同气体模型与低端点恢复','',
  '固定FeII2344、2374、2382，窗口−310到+110km/s，504个全局唯一原生有效像素，expected-fluctuation对角权重。未按拟合残差剪除像素。先运行6、10、14、18、22分量；因为H0AICc持续下降且χ²/ν≈1.61未过门槛，在任何相对位移拟合之前预先限制新增26、30分量、每个三个起点。','',
  '初次扩展错误地让“成功停止”优先于较低的未完成普通端点。这一数值选择问题在输出结论前被独立发现：未完成端点χ²422.187928远低于成功H0χ²463.718227。旧配对和剖面完整保留为before_lower_endpoint_recovery文件；从低端点限定续算，再用相同架构与像素比较H0/H1。','',
  f'最终记录H0χ²={r["null_chi2"]:.9f}/{r["null_nominal_ndf"]}；额外两个相对位移Δχ²={r["delta_chi2"]:.9f}。2382相对2374的D={r["D_m_s"]:+.6f}m/s，局部条件σ={r["conditional_sigma_m_s"]:.6f}m/s。','',
  '共同2382−2374位移接近零；自由模型中的2344相对2374位移约−355.678m/s。因此，两个额外参数合计的Δχ²改善不能直接等同于共同谱线对D的检测。','',
  '这些数值只描述所记录的30分量模型、对角误差、边界和有限优化路径。没有经过完整标定，也没有架构稳定性或全局极小认证。局部σ忽略活动边界的截断和普通气体分量数的不确定性；稀疏profile是真正的干扰参数重优化诊断，不自动构成严格置信区间。五个剖面点中的+2σ点达到300次评价预算而未正常停止，已在图中用叉号标出；中心点目标函数比自由端点低0.000139954，保留为有限优化精度的提示，不按严格最小值或完整置信区间解释。','',
  '独立主复核通过6,469项检查，覆盖本轮四目标全部108次已保存拟合；另有1,247项独立平稳性／局部投影检查。逐目标检查原始FITS、唯一像素、期望涨落误差、CGS单位Voigt模型、卷积、各参数族Jacobian与保存目标。最大导数相对L2误差2.17×10⁻⁸，选择端点及失败诊断的OS21到OS41固定参数χ²变化不超过1.64×10⁻⁵。实现通过复核不代表天体物理模型充分或全局最优。','',
  '文件：`selected_pair.json`、`pair_profile.json`、`lower_endpoint_recovery_protocol.json`、`recovery_outcome.json`、`quadrature_check.json`、`physical_profiles.pdf/.png`均在`results/archive_complete/J225719-100104/`；完整起点、参数、边界、残差、Jacobian、误差和源像素保存在各次JSON/NPZ中。','']
 (ROOT/'reports/archive_complete_j2257_results_cn.md').write_text('\n'.join(targetlines))
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
