#!/usr/bin/env python3
"""Expected-fluctuation-row sensitivity with frozen selected architectures/masks.
The row1 analysis remains preserved. No component-count or window re-selection.
"""
from pathlib import Path
import json
import numpy as np
import archive_profile_pilot as pilot
BASEOUT=pilot.OUT
OUT=BASEOUT/'expected_noise'
BaseModel=pilot.Model
class ExpectedModel(BaseModel):
 def __init__(self,*args,**kwargs):
  super().__init__(*args,**kwargs)
  coadd=np.load(pilot.DATA/f'{self.target}.npz')
  for L in self.lines:
   expected=coadd['expected_fluctuation'][L['source_indices']]
   if not np.all(np.isfinite(expected[L['good']])&(expected[L['good']]>0)):
    raise ValueError('Expected fluctuation missing/nonpositive on frozen retained pixels; cannot compare same mask')
   L['statistical_error']=L['error'].copy()
   L['error']=expected.copy()
  self.cache=None

def run():
 OUT.mkdir(parents=True,exist_ok=True)
 # pilot.fit is reused unchanged; instantiate sensitivity subclass and redirect
 # only its output location. Original Model source and row1 products unchanged.
 pilot.Model=ExpectedModel;pilot.OUT=OUT
 def one(target,n,free,start,suffix):
  r=pilot.fit(target,n,free,start['seed'],start,suffix,max_nfev=600,oversample=33)
  r.update(error_kind='FITS primary array row2: normalized expected fluctuation',sensitivity_script_sha256=pilot.sha(__file__),base_model_script_sha256=pilot.sha(Path(pilot.__file__)),coadd_sha256=pilot.sha(pilot.DATA/f'{target}.npz'),selection_rule='Frozen row1 selected architecture, line windows and valid source pixels; no component-count reselection')
  (OUT/f'{r["name"]}.json').write_text(json.dumps(r,indent=2)+'\n')
  return r
 rows=[]
 for target in pilot.CONFIG:
  selected=json.loads((BASEOUT/f'{target}_selected_refined.json').read_text())
  a=json.loads((BASEOUT/f'{selected["null"]}.json').read_text());b=json.loads((BASEOUT/f'{selected["alternative"]}.json').read_text());n=selected['ncomp']
  h0=one(target,n,False,a,'_expected_os33')
  h1=one(target,n,True,b,'_expected_os33')
  hc=one(target,n,False,h1,'_expected_os33_cross')
  if hc['chi2']<h0['chi2']-1e-4:
   h0=hc
   ac=one(target,n,True,h0,'_expected_os33_cross')
   h1=min([h1,ac],key=lambda x:x['chi2'])
  else:h0=min([h0,hc],key=lambda x:x['chi2'])
  m=ExpectedModel(target,n,False,33)
  row=dict(target=target,ncomp=n,ndata=h0['ndata'],nominal_ndf_null=h0['nominal_ndf'],extra_shifts=h1['npar']-h0['npar'],row1_null=selected['null'],row1_alternative=selected['alternative'],row1_chi2_null=selected['chi2_null'],row1_chi2_alternative=selected['chi2_alternative'],row1_delta_chi2=selected['delta_chi2'],expected_null=h0['name'],expected_alternative=h1['name'],expected_chi2_null=h0['chi2'],expected_chi2_alternative=h1['chi2'],expected_delta_chi2=h0['chi2']-h1['chi2'],both_optimizer_success=h0['optimizer_success'] and h1['optimizer_success'],null_active_bounds=h0['active_bounds'],alternative_active_bounds=h1['active_bounds'],per_line_error_ratio=[dict(line=L['key'],n=int(L['good'].sum()),median_expected_over_statistical=float(np.median(L['error'][L['good']]/L['statistical_error'][L['good']])),min_expected_over_statistical=float(np.min(L['error'][L['good']]/L['statistical_error'][L['good']])),max_expected_over_statistical=float(np.max(L['error'][L['good']]/L['statistical_error'][L['good']]))) for L in m.lines])
  rows.append(row)
 summary=dict(error_kind='FITS primary row2 expected fluctuation',scope='Sensitivity using frozen row1 selected components, lines, windows and masks; not a new model selection and not a covariance model.',comparisons=rows)
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 report=['# 新档案拟合：expected-fluctuation误差行敏感性','',
 'SQUAD合并谱的统计误差（FITS主数组零起始row1）与expected-fluctuation误差（row2）不是同一数组。原档案探索拟合使用统计误差；本附加控制改用expected fluctuation，保留同样气体分量数、谱线、原生像素、窗口与有效掩码，以检查χ²结论对误差选择的依赖。它不替换原结果，不重新选择分量数，也不代表已经处理像素协方差。','',
 '每个目标在33细分像素积分下从已保存的H0/H1参数重新优化，并以H1参数再启动H0。所有原先保留像素上的expected fluctuation均为有限正值，因此无需删点。完整初值、参数、源文件和代码哈希、优化状态保存在本目录JSON/NPZ。','',
 '|目标|相同像素|H0名义自由度|row1 H0 χ²|row1 H1 χ²|row1 Δχ²|row2 H0 χ²|row2 H1 χ²|row2 Δχ²|新增参数|','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
 for r in rows:report.append(f'|{r["target"]}|{r["ndata"]}|{r["nominal_ndf_null"]}|{r["row1_chi2_null"]:.6f}|{r["row1_chi2_alternative"]:.6f}|{r["row1_delta_chi2"]:.6f}|{r["expected_chi2_null"]:.6f}|{r["expected_chi2_alternative"]:.6f}|{r["expected_delta_chi2"]:.6f}|{r["extra_shifts"]}|')
 report += ['', '必须分别解释绝对拟合残差和额外位移改善。某个误差行下的更大χ²不自动证明参数变化；需要继续检查气体结构、仪器响应、连续谱、零点、误差幅度及协方差。若常规模型失配而新增位移也不能充分解决，不能把残差当作α变化或1/ln²时间规律。', '',
 '来源：本项目固定的SQUAD DR1原始FITS；[SQUAD数据论文](https://doi.org/10.1093/mnras/sty2834)。逐行误差比、边界状态、停止条件见`results/archive_expansion/profile_pilot/expected_noise/summary.json`。独立复算结果见同目录`independent_validation.json`。','']
 (pilot.ROOT/'reports/archive_profile_expected_noise_cn.md').write_text('\n'.join(report))
 print(json.dumps(summary,indent=2))
if __name__=='__main__':run()
