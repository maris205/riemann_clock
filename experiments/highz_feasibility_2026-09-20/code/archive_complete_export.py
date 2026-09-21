#!/usr/bin/env python3
"""Stable six-target export; retain all failed-primary diagnostics explicitly."""
import csv,json
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from archive_complete_model import ROOT,OUT,read_json,sha

def per_line_from_npz(path,keys):
 d=np.load(path);rows=[]
 for k in keys:
  g=d[f'{k}_good'];r=(d[f'{k}_model'][g]-d[f'{k}_flux'][g])/d[f'{k}_error'][g]
  rows.append({'line':k,'pixels':int(g.sum()),'chi2':float(r@r),'chi2_per_pixel':float(r@r/g.sum())})
 return rows

def main():
 recovery=OUT/'J225719-100104/recovery_outcome.json'
 if not recovery.exists() or read_json(recovery).get('status')!='recovered_pair_saved':
  raise RuntimeError('J225719 recovery not finalized; refusing to export obsolete pair')
 commonpath=ROOT/'results/common_pair/refined_summary.json';common=read_json(commonpath);rows=[]
 for a in common['comparisons']:
  f=a['full_comparison'];nu=a['ndata']-f['null_npar'];h0path=ROOT/f['null_file'];h0=read_json(h0path)
  line=per_line_from_npz(h0path.with_suffix('.npz'),a['lines'])
  gate=f['chi2_null']/nu<=1.5 and max(l['chi2_per_pixel'] for l in line)<=1.8
  rows.append({'target':a['target'],'z_abs':a['z_abs'],'pixels':a['ndata'],'ncomp':a['ncomp'],
   'Q0':f['chi2_null'],'Q1':f['chi2_alternative'],'nu0':nu,'nu1':a['ndata']-f['alternative_npar'],
   'delta_Q':f['delta_chi2'],'extra_shifts':f['added_parameters'],
   'diagnostic_H1_status':'conditional_common_pair_companion_complete',
   'D_m_s':a['D_m_s'],'conditional_local_error_m_s':a['conditional_sigma_m_s'],
   'diagnostic_only_D_m_s':None,'adequacy_flag':gate,
   'adequacy_reasons':['Retrospective same numeric residual screen for earlier compact pilots; no claim of prospective six-target gate.',
      'Fixed architecture and bounds; derivative error and bounded profile differ in nonlinear/bound-active cases.'],
   'null_per_line':line,'both_optimizer_stopped':f['null_optimizer_success'] and f['alternative_optimizer_success'],
   'source_summary':str(commonpath.relative_to(ROOT)),'source_summary_sha256':sha(commonpath),
   'null_fit_file':f['null_file'],'alternative_fit_file':f['alternative_file'],
   'profile_intervals':a['profile_intervals'],'physical_precision_measurement':False,
   'first_order_stationarity_certified':False})
 for target in ['J053007-250329','J225719-100104']:
  folder=OUT/target;p=read_json(folder/'selected_pair.json');null=read_json(folder/(p['null']+'.json'));alt=read_json(folder/(p['alternative']+'.json'))
  local=p['local_pair_uncertainty'];profile=read_json(folder/'pair_profile.json')
  rows.append({'target':target,'z_abs':p['z_abs'],'pixels':p['ndata'],'ncomp':p['ncomp'],
   'Q0':p['chi2_null'],'Q1':p['chi2_alternative'],'nu0':p['nominal_ndf_null'],'nu1':alt['nominal_ndf'],
   'delta_Q':p['delta_chi2'],'extra_shifts':p['extra_shifts'],
   'diagnostic_H1_status':('conditional_H1_complete_sparse_profile_partly_capped' if any(not q['optimizer_success'] for q in profile['grid']) else 'conditional_H1_and_sparse_profile_completed'),
   'D_m_s':local['shift_m_s'],'conditional_local_error_m_s':local['rank_sensitivity'][1]['conditional_sigma_m_s'],
   'diagnostic_only_D_m_s':None,'adequacy_flag':p['ordinary_gate'],
   'adequacy_reasons':['Passed prespecified residual screen under expected diagonal errors; not physical or population-null certification.',
     'Gas count is final search upper edge; finite optimizer paths, bound-active nuisance and calibration remain unresolved.'],
   'null_per_line':null['per_line'],'both_optimizer_stopped':p['both_optimizer_success'],
   'source_summary':str((folder/'selected_pair.json').relative_to(ROOT)),'source_summary_sha256':sha(folder/'selected_pair.json'),
   'null_fit_file':str((folder/(p['null']+'.json')).relative_to(ROOT)),
   'alternative_fit_file':str((folder/(p['alternative']+'.json')).relative_to(ROOT)),
   'sparse_profile':profile['grid'],
   'sparse_profile_all_optimizer_stopped':all(q['optimizer_success'] for q in profile['grid']),
   'sparse_profile_lowest_minus_free_Q':min(q['delta_from_free'] for q in profile['grid']),
   'sparse_profile_is_confidence_interval':False,'physical_precision_measurement':False,
   'first_order_stationarity_certified':False})
 # Failed primary targets remain in the same sample and Q table.
 folder=OUT/'J064326-504112';p=read_json(folder/'failure_pair/summary.json')
 rows.append({'target':p['target'],'z_abs':2.659,'pixels':p['ndata'],'ncomp':p['ncomp'],
  'Q0':p['chi2_null'],'Q1':p['chi2_alternative'],'nu0':p['nominal_ndf_null'],'nu1':p['nominal_ndf_alternative'],
  'delta_Q':p['delta_chi2'],'extra_shifts':p['extra_shift_parameters'],
  'diagnostic_H1_status':'failed_primary_model_H1_did_not_repair_weak_1611_mismatch',
  'D_m_s':None,'conditional_local_error_m_s':None,'diagnostic_only_D_m_s':p['shifts_m_s'].get('2382'),
  'adequacy_flag':False,'adequacy_reasons':['Four-line primary FeII1611 chi2/pixel exceeds1.8; remainsabove1.8underH1.',
      'Weak-line quality/gas/continuum ambiguity not isolated; sky flag alone does not explain largest residual.'],
  'null_per_line':p['null_per_line'],'both_optimizer_stopped':p['both_optimizer_success'],
  'source_summary':str((folder/'failure_pair/summary.json').relative_to(ROOT)),'source_summary_sha256':sha(folder/'failure_pair/summary.json'),
  'null_fit_file':str((folder/'failure_pair'/p['null_relative_path']).resolve().relative_to(ROOT)),
  'alternative_fit_file':str((folder/'failure_pair'/p['alternative_relative_path']).resolve().relative_to(ROOT)),
  'physical_precision_measurement':False,'first_order_stationarity_certified':False})
 folder=OUT/'J233156-090802';p=read_json(folder/'failure_pair/summary.json');n=p['null'];a=p['alternative']
 rows.append({'target':p['target'],'z_abs':2.143,'pixels':p['native_pixels'],'ncomp':p['ncomp'],
  'Q0':n['chi2'],'Q1':a['chi2'],'nu0':n['nominal_ndf'],'nu1':a['nominal_ndf'],
  'delta_Q':p['delta_chi2'],'extra_shifts':p['extra_region_shifts'],
  'diagnostic_H1_status':'failed_primary_model_H1_did_not_repair_total_and_2382_mismatch',
  'D_m_s':None,'conditional_local_error_m_s':None,'diagnostic_only_D_m_s':a['shifts_m_s'].get('2382'),
  'adequacy_flag':False,'adequacy_reasons':['Total chi2/nu>1.5 and2382 chi2/pixel>1.8 inH0andH1.',
      'Overlapping1611 opacity modeledinside1608region without duplicatepixels; remaining gas/model mismatch unresolved.'],
  'null_per_line':n['per_line'],'both_optimizer_stopped':p['both_selected_optimizer_success'],
  'source_summary':str((folder/'failure_pair/summary.json').relative_to(ROOT)),'source_summary_sha256':sha(folder/'failure_pair/summary.json'),
  'null_fit_file':str((folder/(n['name']+'.json')).relative_to(ROOT)),
  'alternative_fit_file':str((folder/'failure_pair/fits'/p['target']/(a['name']+'.json')).relative_to(ROOT)),
  'physical_precision_measurement':False,'first_order_stationarity_certified':False})
 rows.sort(key=lambda r:r['z_abs'])
 for row in rows:
  for field in ['null_fit_file','alternative_fit_file']:
   if not (ROOT/row[field]).is_file():raise FileNotFoundError(row[field])
   row[field+'_sha256']=sha(ROOT/row[field])
 result={'created_utc':datetime.now(timezone.utc).isoformat(),
   'quantity':'D=velocity(FeII2382)-velocity(FeII2374), conditional relative registration; not alpha.',
   'error_model':'All Q comparisons use native expected-fluctuation arrays as diagonalweights.',
   'selection_warning':'All six automatic candidates retained. Adequacygate can reject realfrequencychanges; passedsubsetcannot establish populationnull. Failurepairshift values are diagnostic-only, with mainD/errors null.',
   'calibration_warning':'No row is a calibrated cosmological precision measurement; first-order/globalstationarity is not certified.',
   'rows':rows,'generator_sha256':sha(__file__)}
 (OUT/'six_target_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
 keys=['target','z_abs','pixels','ncomp','Q0','Q1','nu0','nu1','delta_Q','extra_shifts','diagnostic_H1_status','D_m_s','conditional_local_error_m_s','diagnostic_only_D_m_s','adequacy_flag','both_optimizer_stopped','source_summary']
 with (OUT/'six_target_comparison.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows([{k:r[k] for k in keys} for r in rows])
 report=['# 六个自动候选的完整状态表','',
  '这六个对象均完成了固定主像素和 expected-fluctuation 对角误差下的共同气体模型 H0 与允许逐谱区相对位移的 H1 比较。四个对象通过数值残差门槛，可提供依赖模型的条件位移描述；另外两个普通模型仍不充分，保留在表中，仅作失配诊断。任何对象都没有被认定为已校准的宇宙学精密测量。','',
  '|目标|吸收红移|原生像素|气体分量|Q0 / ν0|Q1 / ν1|ΔQ / 额外位移|状态|',
  '|---|---:|---:|---:|---:|---:|---:|---|']
 for r in rows:
  report.append(f'|{r["target"]}|{r["z_abs"]:.3f}|{r["pixels"]}|{r["ncomp"]}|{r["Q0"]:.6f} / {r["nu0"]}|{r["Q1"]:.6f} / {r["nu1"]}|{r["delta_Q"]:.6f} / {r["extra_shifts"]}|'+('条件位移描述' if r['adequacy_flag'] else '失配诊断；不报告精密位移')+'|')
 report += ['', '共同描述量为 D=v(Fe II 2382)−v(Fe II 2374)。下列局部误差已投影掉其余全部拟合干扰方向，但不包含波长标定、气体架构选择、完整协方差或活动边界效应。','',
  '|目标|D (m/s)|局部条件 σ (m/s)|','|---|---:|---:|']
 for r in rows:
  if r['D_m_s'] is not None:report.append(f'|{r["target"]}|{r["D_m_s"]:+.3f}|{r["conditional_local_error_m_s"]:.3f}|')
 report += ['',
  'J2257 的合计 ΔQ=13.53 对应两个位移参数，不能直接解释成共同谱线对 D 的检测；其 2344 相对位移约−356m/s，而 D 接近零。它的五点稀疏剖面中 +2σ 点达到评价预算，中心点比自由端点低约1.40×10⁻⁴；剖面保留为数值诊断，不给出严格置信区间。','',
  'J0643 的 1611 谱线失配和 J2331 的总体／2382 谱线失配，在加入自由位移后仍未修复。质量门槛也可能排除真实频率变化，所以通过门槛的子样本不能建立总体不变结论。四个新增目标均达到此次气体分量数搜索上限，且所有目标仍有标定和局部最优限制；不能把这些位移直接转换成 α 或 1/ln²(t/t*) 的物理约束。','',
  '机器可读结果：`results/archive_complete/six_target_comparison.json` 和 `.csv`；JSON 包含每个源文件路径、SHA-256、逐线残差、优化状态、主分析与失败诊断的区分。旧的两个紧凑对象使用 `results/common_pair/refined_summary.json` 中重新剖面后的最新有限端点，未覆盖旧历史结果。','']
 (ROOT/'reports/archive_complete_six_target_summary_cn.md').write_text('\n'.join(report))
 print(json.dumps([{k:r[k] for k in ['target','z_abs','Q0','Q1','delta_Q','extra_shifts','D_m_s','conditional_local_error_m_s','adequacy_flag']} for r in rows],indent=2))
if __name__=='__main__':main()
