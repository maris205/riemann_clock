#!/usr/bin/env python3
"""Bounded free-shift diagnostic for a primary ordinary-model gate failure.

No precision descriptor is promoted. Same fixed three regions, 1176 unique
native pixels, 26 gas components, expected diagonal errors, Fe1611 overlap.
The primary gate is a quality screen, not a population-wide null test.
"""
from pathlib import Path
from datetime import datetime,timezone
import json
import numpy as np
from archive_complete_j2331_model import NeighborModel,engine
ROOT=engine.ROOT
PRIMARY=ROOT/'results/archive_complete/J233156-090802'
OUT=PRIMARY/'failure_pair'
def dump(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 configfile=PRIMARY/'selection.json';config=engine.read_json(configfile)
 ordinary=engine.read_json(PRIMARY/'ordinary_model_summary.json')
 startfile=PRIMARY/(ordinary['selected_null']+'.json');start=engine.read_json(startfile)
 assert config['neighbor_opacity_adapter_sha256']==engine.sha(ROOT/config['neighbor_opacity_adapter'])
 assert start['ncomp']==26 and start['ndata']==1176 and not ordinary['passes_conditional_adequacy_gate']
 m0=NeighborModel(config,26,False,21);m1=NeighborModel(config,26,True,21)
 ids0=np.concatenate([L['source_indices'][L['good']] for L in m0.lines]);ids1=np.concatenate([L['source_indices'][L['good']] for L in m1.lines])
 assert np.array_equal(ids0,ids1) and len(np.unique(ids0))==len(ids0)==1176
 protocol={'frozen_utc':datetime.now(timezone.utc).isoformat(),'authorization':'Bounded failure-case diagnostic requested by root after primary adequacy gate failed; amendment recorded before additional fits.',
 'scientific_reason':'A conventional-model quality gate must not prevent testing whether free shifts can repair its gross mismatch, or be mistaken for a population-wide no-variation conclusion.',
 'primary_disposition_preserved':'Primary descriptor D remains not estimated; any diagnostic shift stays nonprecision regardless of delta_chi2.',
 'target':config['target'],'configuration_sha256':engine.sha(configfile),'configuration':config,'initial_primary_null_file':str(startfile.relative_to(ROOT)),'initial_primary_null_sha256':engine.sha(startfile),
 'initial_primary_null_npz_sha256':engine.sha(startfile.with_suffix('.npz')),'initial_primary_null_chi2':start['chi2'],'ncomp':26,'native_pixels':1176,'unique_source_indices':len(np.unique(ids0)),
 'same_pixels_H0_H1':True,'neighbor_opacity':'Same-gas1611 included only inside1608 region; no extra likelihood region.',
 'oversample':21,'first_H1_max_nfev':600,'H0_cross_max_nfev':500,'optional_H1_cross_max_nfev':500,
 'optional_cross_trigger':'Run one H1 cross if H0 cross improves primary H0 by more than 0.001 in chi2.',
 'profile_uncertainty':False,'free_shift_bounds_km_s':[-1.,1.],'model_adapter_sha256':engine.sha(ROOT/'code/archive_complete_j2331_model.py'),'shared_engine_sha256':engine.sha(ROOT/'code/archive_complete_model.py'),'runner_sha256':engine.sha(__file__),
 'limitations':['Exploratory and prompted by observed ordinary-model inadequacy.','No calibrated discovery significance, precision line-shift datum, alpha inference or cosmic-age-law fit.','Calibration, blends, gas architecture, active boundaries and local minima remain unresolved.']}
 protocolpath=OUT/'protocol.json'
 if protocolpath.exists():raise RuntimeError('Protocol already exists: do not overwrite a completed or partial diagnostic.')
 dump(protocolpath,protocol)
 engine.Model=NeighborModel;engine.OUT=OUT/'fits'
 alt=engine.fit(config,26,True,start['seed'],start,name='diagnostic_H1_from_primary_os21',oversample=21,max_nfev=600)
 cross=engine.fit(config,26,False,start['seed'],alt,name='diagnostic_H0_cross_os21',oversample=21,max_nfev=500)
 attempts=[alt,cross];null=min([start,cross],key=lambda r:r['chi2'])
 if cross['chi2']<start['chi2']-.001:
  alt2=engine.fit(config,26,True,start['seed'],cross,name='diagnostic_H1_cross_os21',oversample=21,max_nfev=500);attempts.append(alt2)
  alt=min([alt,alt2],key=lambda r:r['chi2'])
 def compact(r):
  return {k:r[k] for k in ['name','chi2','chi2_per_ndf','ndata','npar','nominal_ndf','optimizer_success','optimizer_status','nfev','stationarity','active_bounds','shifts_m_s','per_line']}
 summary={'target':config['target'],'primary_adequacy_gate':False,'primary_precision_descriptor_D':'not estimated; unchanged',
 'diagnostic_only':True,'ncomp':26,'native_pixels':1176,'null':compact(null),'alternative':compact(alt),'attempts':[compact(r) for r in attempts],
 'delta_chi2':null['chi2']-alt['chi2'],'extra_region_shifts':alt['npar']-null['npar'],'both_selected_optimizer_success':bool(null['optimizer_success'] and alt['optimizer_success']),
 'alternative_meets_same_numeric_quality_thresholds':bool(alt['optimizer_success'] and alt['chi2_per_ndf']<=1.5 and max(r['chi2_per_pixel'] for r in alt['per_line'])<=1.8),
 'interpretation':'Bounded gross-mismatch repair diagnostic only. Original gate is a quality screen, not a population-wide absence-of-variation test. Diagnostic offsets must not be promoted to precision D regardless of improvement.',
 'protocol_sha256':engine.sha(protocolpath),'primary_configuration_sha256':engine.sha(configfile)}
 dump(OUT/'summary.json',summary)
 lines=['# J233156−090802：常规模型门槛未通过后的有限自由位移诊断','',
 '原有门槛筛选保留为质量判断；它不能排除真实频移本身造成常规模型失配，也不能被解释为整个候选样本均无变化。因此在追加拟合前冻结此项检验，只判断两个额外谱区位移能否改善原有明显失配。','',
 '仍使用原三谱区、1176个全局唯一源像素、26个共享气体分量、expected-fluctuation对角误差、OS21及1608谱区内同气体1611不透明度。没有剪除像素、改变线集或重新选择架构。','',
 '| 条件比较 | χ² | 名义自由度 | χ²/ν | 成功停止 |','|---|---:|---:|---:|---|',
 f"| 保留的H0 | {null['chi2']:.9f} | {null['nominal_ndf']} | {null['chi2_per_ndf']:.6f} | {null['optimizer_success']} |",
 f"| 诊断H1 | {alt['chi2']:.9f} | {alt['nominal_ndf']} | {alt['chi2_per_ndf']:.6f} | {alt['optimizer_success']} |",'',
 f"增加{summary['extra_region_shifts']}个谱区位移的条件改善为Δχ²={summary['delta_chi2']:.9f}。H1是否达到同一组数值质量阈值：{summary['alternative_meets_same_numeric_quality_thresholds']}。",'',
 '**原主分析中的精密D测量仍为未估计。** 无论上述改善多少，诊断位移不进入宇宙年龄曲线、不转换成α，也不附加高斯发现显著性。1608位移作用于包含1608和1611的整个局部谱区；目标2382−2374仍只作为非精密诊断。', '',
 f"保留H0活动边界：`{null['active_bounds']}`；H1活动边界：`{alt['active_bounds']}`。一阶指标与每次停止状态完整保存在summary及各拟合JSON。成功ftol/xtol停止不保证全局最优。",'',
 '该检验不重新定义原质量门槛。气体结构、局部极小值、仪器响应、波长标定、误差协方差和额外混叠仍需独立处理。']
 (ROOT/'reports/archive_complete_j2331_failure_pair_cn.md').write_text('\n'.join(lines)+'\n')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
