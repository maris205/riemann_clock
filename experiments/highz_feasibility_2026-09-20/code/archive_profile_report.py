#!/usr/bin/env python3
"""Export completed archival physical-profile pilot without rerunning fits."""
import json
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from archive_profile_pilot import OUT,ROOT,CONFIG
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
selected=[]
with PdfPages(OUT/'archive_profile_fits.pdf') as pdf:
 for target,cfg in CONFIG.items():
  s=json.loads((OUT/f'{target}_selected_refined.json').read_text());selected.append(s)
  h0=json.loads((OUT/f'{s["null"]}.json').read_text());h1=json.loads((OUT/f'{s["alternative"]}.json').read_text())
  d=np.load(OUT/f'{s["null"]}.npz');alt=np.load(OUT/f'{s["alternative"]}.npz');keys=cfg['lines']
  fig,axs=plt.subplots(len(keys),2,figsize=(10.5,2.0*len(keys)),squeeze=False,gridspec_kw={'width_ratios':[1.4,1]})
  for i,k in enumerate(keys):
   v=d[f'{k}_v'];f=d[f'{k}_flux'];e=d[f'{k}_error'];m=d[f'{k}_model'];g=d[f'{k}_good']
   axs[i,0].errorbar(v[g],f[g],e[g],fmt='.',ms=3,color='#253c50',alpha=.8,label='Native coadd pixels')
   axs[i,0].plot(v,m,c='#007f79',lw=1.5,label='Shared gas; fixed line frequencies')
   axs[i,0].plot(v,alt[f'{k}_model'],c='#d77622',lw=1,ls='--',label='With relative line shifts')
   axs[i,0].set_ylabel(f'Fe II {k}\nNormalized flux');axs[i,0].set_ylim(-.06,1.09)
   axs[i,1].plot(v[g],(f[g]-m[g])/e[g],'.-',lw=.6,ms=3,color='#253c50');axs[i,1].axhline(0,c='k',lw=.6)
   axs[i,1].axhspan(-2,2,color='#007f79',alpha=.09);axs[i,1].set_ylabel('Null residual / error');axs[i,1].set_ylim(-5,5)
  axs[0,0].legend(fontsize=7,loc='lower left');axs[-1,0].set_xlabel('Velocity relative to catalogue redshift (km/s)');axs[-1,1].set_xlabel('Velocity (km/s)')
  fig.suptitle(f'{target}, z = {cfg["z"]:.3f}; {s["ncomp"]} shared gas components\nNull chi-square = {s["chi2_null"]:.2f} / {s["nominal_ndf_null"]} nominal dof; extra-shift improvement = {s["delta_chi2"]:.2f} ({s["extra_shifts"]} parameters)\nExploratory profile fit; diagonal pixel errors; no alpha or time-law measurement',fontsize=11)
  fig.tight_layout(rect=(0,0,1,.94));pdf.savefig(fig);fig.savefig(OUT/f'{target}_physical_profiles.png',dpi=170);fig.savefig(OUT/f'{target}_physical_profiles.pdf');plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(10,6.6))
for j,(target,cfg) in enumerate(CONFIG.items()):
 rows=json.loads((OUT/f'{target}_comparison.json').read_text())['comparisons'];ns=[r['ncomp'] for r in rows]
 red=[r['chi2_null']/r['nominal_ndf_null'] for r in rows];deltas=[r['delta_chi2'] for r in rows]
 axs[0,j].plot(ns,red,'o-',c='#007f79');axs[0,j].axhline(1,c='gray',ls='--',lw=.8);axs[0,j].set_yscale('log');axs[0,j].set_title(f'{target}\nz = {cfg["z"]:.3f}');axs[0,j].set_ylabel('Null chi-square / nominal dof')
 axs[1,j].plot(ns,deltas,'o-',c='#d77622');axs[1,j].set_yscale('log');axs[1,j].set_ylabel('Extra-shift chi-square improvement');axs[1,j].set_xlabel('Number of shared gas components')
 for a in axs[:,j]:a.set_xticks(ns);a.grid(alpha=.15)
fig.suptitle('Ordinary gas structure must be fitted before interpreting relative shifts\nExploratory same-pixel comparisons; diagonal errors; no detection significance',fontsize=12);fig.tight_layout(rect=(0,0,1,.92));fig.savefig(OUT/'archive_profile_summary.png',dpi=180);fig.savefig(OUT/'archive_profile_summary.pdf');plt.close(fig)
lines=['# 两个新增档案吸收系统：实际共享气体谱线拟合','',
'本次把档案候选推进到实际物理线型拟合。两个目标都能在探索性的常规共享气体模型下得到可接受的残差；额外自由相对位移的改善有限。这是分析可行性和常规模型核查的结果，不是精细结构常数变化或对数时间规律的测量。','',
'## 固定分析范围与物理模型','',
'数据来自本项目已归档的 [UVES SQUAD DR1](https://doi.org/10.1093/mnras/sty2834) 合并谱及原生像素窗口。目标和窗口在查看谱形之后选择，属于探索性分析。J232128−105122（z=1.629）使用 Fe II 2260、2344、2374、2382、2586，固定速度窗口 −30 到 +35 km/s；窗口避开 2344 更红侧的额外结构。J004131−493611（z=2.248）使用 Fe II 1608、2374、2382，固定窗口 −30 到 +85 km/s。速度零点由目录中四舍五入的吸收红移定义，实际共同速度通过气体分量自由拟合。没有依据拟合残差删除像素。','',
'每个气体分量具有共享的 Fe II 柱密度 N、速度 v 和 Doppler 参数 b。原子同位素波长、自然展宽及已经按丰度加权的振子强度与 ESPRESSO 检验采用同一原子文件，原子数据背景见 [Murphy & Berengut (2014)](https://doi.org/10.1093/mnras/stt2204)。先对总光学深度取 exp(−τ)，再卷积 Gaussian 仪器响应并进行像素积分；实际网格步长直接来自波长数组，不用四舍五入的 UP_DISP。没有把 ESPRESSO 作者的零点或连续谱修正搬到 UVES。','',
'两种假设都拟合每条线的连续谱常数/斜率、加性零点及仪器 FWHM。零点限制±0.02，连续谱常数±0.05、每100 km/s斜率±0.10，FWHM为目录名义值的0.6–1.4倍。替代模型仅增加各线相对2374的位移，范围±1 km/s。两假设使用完全相同的像素和误差；误差暂按对角处理，没有测量这些新目标的像素协方差。','',
'运行1、2、3、4、6个共享气体分量，每种至少三个初始位置/柱密度扰动，并进行H0/H1交叉初始值拟合。对J004131额外增加普通共享弱分量的初始值于−10 km/s；原先所有通用初始值都集中在主吸收区，未找到这个常规模型极小值。所有尝试都保留。按H0的探索性对角误差AICc选择展示模型后，像素积分从11细分提高到33细分，并重新优化和交叉拟合；这不是盲化预注册的模型选择。','',
'## 完成的拟合结果','',
'| 目标 | 分量 | 像素 | 常规模型 χ² / 名义自由度 | 加相对位移后的 χ² | 新参数 | Δχ² |','|---|---:|---:|---:|---:|---:|---:|']
for s in selected:
 lines.append(f'| {s["target"]} | {s["ncomp"]} | {s["ndata"]} | {s["chi2_null"]:.6f} / {s["nominal_ndf_null"]} | {s["chi2_alternative"]:.6f} | {s["extra_shifts"]} | {s["delta_chi2"]:.6f} |')
lines += ['',
'J232128的常规模型本身已能解释所选线型；3/4/6分量时，额外位移的Δχ²约6.84、6.61、5.50。弱2260线的位移并不精确，不应因为点估计较大而声称异常。', '',
'J004131最初的6分量常规模型χ²约182.31，其中2382在−20到0 km/s有明显剩余吸收。将一个普通、与其余跃迁共享的弱气体分量初始位置放到这一范围，再优化后χ²降到113.19。没有添加只作用于2382的异常频率机制，也没有剪掉该区域。这是为什么先检查气体结构与多初值收敛很重要的一个具体例子。', '',
'两目标在最终展示模型中仍有加性零点触及允许边界（J232128的2260、2374；J004131的2374）。这说明零点/连续谱/气体参数仍有退化，优化停止成功不能当作所有系统误差已经确定。自由位移结果保存在JSON用于可复算，但本报告不把它们当作校准完成的速度测量，也不转换成α变化。','',
'## 分量数敏感性（11细分初始探索；最终展示值采用33细分）','',
'| 目标 | 分量数 | H0 χ²/名义自由度 | Δχ² | 额外位移参数数 |','|---|---:|---:|---:|---:|']
for target,cfg in CONFIG.items():
 for r in json.loads((OUT/f'{target}_comparison.json').read_text())['comparisons']:
  lines.append(f'| {target} | {r["ncomp"]} | {r["chi2_null"]:.2f}/{r["nominal_ndf_null"]} | {r["delta_chi2"]:.2f} | {len(cfg["lines"])-1} |')
lines += ['',
'少分量模型尤其是J004131单分量，给出很大的位移改善，但绝对线型拟合严重失败，位移还碰到边界。这些结果被保留为模型失配示例，不能用作变化信号。增加普通气体结构后，大幅改善消失。','',
'## 结论与仍需处理的项目','',
'新增档案数据确实能够做实际共享物理线型拟合，但本次没有提供比原HE0515结果更强的变化证据。两个新增系统都支持先用常规气体、仪器及误差模型解释数据的工作顺序。其余档案候选仍需按质量筛查逐一推进，不应把目录中的目标数直接计为有效时间点。','',
'还需要每个目标自己的噪声协方差、曝光层面的标定和合并核查、连续谱/零点/LSF敏感性、更完整的分量搜索和邻近混合线识别。多初值只提高可靠性，不证明全局最优。不同红移的这两个探索性模型不能直接拼成1/ln²(t/t*)曲线；必须先明确原子响应并获得可比较的物理观测量。','',
'## 可复算文件','',
'- `code/archive_profile_pilot.py`：物理模型和初始多分量拟合。',
'- `code/archive_profile_refine.py`：记录新增普通气体初始值、模型选择与33细分重拟合。',
'- `code/archive_profile_report.py`：本报告和图片导出。',
'- `results/archive_expansion/profile_pilot/policy.json`：冻结窗口和参数边界。',
'- `results/archive_expansion/profile_pilot/*_selected_refined.json`：最终选择、χ²、位移及优化状态。',
'- `results/archive_expansion/profile_pilot/*.json/.npz`：全部初始值、参数、残差、Jacobian、源文件和原子文件SHA256。',
'- `results/archive_expansion/profile_pilot/archive_profile_fits.pdf`：两个目标的实际光谱、物理模型、标准化残差。',
'- `results/archive_expansion/profile_pilot/archive_profile_summary.png`：分量数依赖。',
'- `results/archive_expansion/profile_pilot/independent_validation.json` 和 `reports/archive_profile_validation.md`：独立实现核查：67个已保存拟合、3,362项检查全部通过，包含最终33细分选择结果；33到55细分的Δχ²变化低于3×10⁻⁸。','']
expected_path=OUT/'expected_noise/summary.json'
if expected_path.exists():
 expected=json.loads(expected_path.read_text())
 lines += ['## 统计误差与expected-fluctuation误差的追加控制','',
 '上面的主分析使用FITS主数组统计误差row1。SQUAD还提供expected-fluctuation误差row2；我们在相同分量数、谱线、窗口及有效像素下，完成了33细分的成对重拟合及H0交叉初值控制，没有重新选择分量数。全部原选定像素上的row2误差均为有限正值。','',
 '|目标|row1 Δχ²|row2 H0 χ² / 名义自由度|row2 H1 χ²|row2 Δχ²|','|---|---:|---:|---:|---:|']
 for row in expected['comparisons']:
  lines.append(f'|{row["target"]}|{row["row1_delta_chi2"]:.6f}|{row["expected_chi2_null"]:.6f}/{row["nominal_ndf_null"]}|{row["expected_chi2_alternative"]:.6f}|{row["expected_delta_chi2"]:.6f}|')
 lines += ['', '采用expected-fluctuation误差后，两个系统的额外位移改善都略减小，未改变本次的条件结论。两种误差行都没有替代像素协方差或仪器标定；因此仍不报告变化显著性。详细过程见`reports/archive_profile_expected_noise_cn.md`，所有六个追加拟合及独立验证保存在`results/archive_expansion/profile_pilot/expected_noise/`。','']
validation_path=OUT/'expected_noise/independent_validation.json'
if validation_path.exists():
 validation=json.loads(validation_path.read_text())
 if validation.get('passed'):
  lines += ['追加expected-fluctuation控制的独立验证：6个新拟合、615项检查全部通过。验证从原始FITS直接读取row2，复核与主分析相同的像素/架构，并独立复算权重后的残差、Jacobian和χ²；详见`reports/archive_profile_expected_validation.md`。','']
  er=ROOT/'reports/archive_profile_expected_noise_cn.md'
  if er.exists() and '615项检查全部通过' not in er.read_text():
   er.write_text(er.read_text()+'\n独立验证完成：6个新增拟合、615项检查全部通过。从原始FITS直接验证row2误差、保持不变的像素/架构，以及重新加权的残差、Jacobian和χ²。报告见`reports/archive_profile_expected_validation.md`。\n')
(ROOT/'reports/archive_profile_results_cn.md').write_text('\n'.join(lines))
summary=dict(selected=selected,conclusion='Both exploratory ordinary shared-gas models fit their fixed native-pixel windows. Additional relative shifts provide modest conditional improvements only; no alpha/time-law measurement. Missed conventional gas start explained initial J004131 residual.',report='reports/archive_profile_results_cn.md')
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
