#!/usr/bin/env python3
"""Export already-computed conditional common-pair profiles; no fitting."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/common_pair'
summary=json.loads((OUT/'refined_summary.json').read_text());rows=summary['comparisons']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def interval(row,level):return next(x for x in row['profile_intervals'] if abs(x['delta_chi2_level']-level)<1e-10)
def endpoints(c):return [c[k]['crossing_m_s'] if c[k] else None for k in ['lower','upper']]
def fpair(c):
 p=endpoints(c);return '未封闭' if None in p else f'[{p[0]:.1f}, {p[1]:.1f}]'
review_path=OUT/'review/validation.json';review=json.loads(review_path.read_text()) if review_path.exists() else None
harmonized=[]
for row in rows:
 h={k:row[k] for k in ['target','z_abs','D_m_s','conditional_sigma_m_s','profile_grid_D_m_s','profile_delta_chi2','active_bounds','profile_all_selected_terminated','profile_contours_closed','profile_baseline_resolved_to_tolerance','unrestricted_reference_gap_chi2','point_estimate_reference_file','point_estimate_reference_sha256','ndata','ncomp','lines','oversample','sourcefit_paths','sourcefit_hashes','full_comparison']}
 h.update(observable='D = v(Fe II 2382)-v(Fe II 2374)',unit='m/s',anchor_line=2374,target_line=2382,log_energy_ratio_sign='delta ln(E2382/E2374) = -D/c; kinematic ratio only, not an alpha or clock response',local_sigma_role='Unconstrained nuisance-tangent curvature in row2 diagonal errors; not reduced-chi-square rescaled or coverage calibrated',profile_delta1_interval_m_s=endpoints(interval(row,1.)),profile_delta3p84_interval_m_s=endpoints(interval(row,3.84)),profile_is_conditional_connected_contour=True,precision_measurement=False,physical_time_inference=False,other_relative_line_shifts_float=[k for k in row['lines'] if k not in [2374,2382]],error_kind=row['error_kind'],condition='Fixed gas architecture, terrestrial isotope mixture, bounded Gaussian LSF/continuum/zero, diagonal row2 errors; missing calibration/covariance/architecture uncertainty',source_summary_file='results/common_pair/refined_summary.json',source_summary_sha256=sha(OUT/'refined_summary.json'))
 h['profile_crossing_brackets']={str(level):{side:None if interval(row,level)[side] is None else dict(coordinate_bracket_m_s=interval(row,level)[side]['bracket_m_s'],crossing_m_s=interval(row,level)[side]['crossing_m_s'],endpoint_fits=interval(row,level)[side]['fit_files'],endpoint_order='inside-to-outside for fit files; coordinate bracket independently sorted') for side in ['lower','upper']} for level in [1.,3.84]}
 harmonized.append(h)
(OUT/'harmonized_measurements.json').write_text(json.dumps(dict(scope='Conditional common-observable summaries for feasibility/identifiability only; no accepted cosmological precision measurements.',comparisons=harmonized),indent=2)+'\n')
plt.rcParams.update({'font.size':9,'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':7})
fig,axes=plt.subplots(1,2,figsize=(9,4.15),sharey=True)
for ax,row in zip(axes,rows):
 x=np.array(row['profile_grid_D_m_s']);y=np.array(row['profile_delta_chi2']);at=row['D_m_s'];sigma=row['conditional_sigma_m_s']
 order=np.argsort(x);x=x[order];y=y[order]
 ax.plot(x,y,'o-',color='#2166ac',markersize=3.2,lw=1.2,label='Bounded nuisance profile')
 grid=np.linspace(x.min(),x.max(),600)
 ax.plot(grid,((grid-at)/sigma)**2,color='#d38a2f',ls='--',lw=1.2,label='Unconstrained local tangent')
 for level,style in [(1.,':'),(3.84,'-.')]:ax.axhline(level,color='.45',ls=style,lw=.85)
 ax.axvline(0,color='.65',lw=.8);ax.axvline(at,color='#2166ac',lw=.6,alpha=.65)
 lo,hi=endpoints(interval(row,1.))
 if lo is not None and hi is not None:ax.axvspan(lo,hi,color='#2166ac',alpha=.08)
 ax.set_ylim(-.12,10);ax.set_xlim(x.min()-.03*np.ptp(x),x.max()+.03*np.ptp(x))
 ax.set_title(f'{row["target"]}\n$z={row["z_abs"]:.3f}$; {row["ncomp"]} components, {len(row["lines"])} lines')
 ax.set_xlabel(r'$D=v_{2382}-v_{2374}$ (m s$^{-1}$)');ax.grid(alpha=.15)
 ax.legend(loc='upper center',frameon=True)
axes[0].set_ylabel(r'$\Delta\chi^2$ relative to best recorded endpoint')
fig.suptitle('Common Fe II pair: conditional profile under fixed model and error assumptions',fontsize=11)
fig.text(.5,.008,'Other relative line shifts float. Contours at 1 and 3.84 are reference objective levels, not calibrated coverage.',ha='center',fontsize=7.5)
fig.tight_layout(rect=[0,.035,1,.94]);fig.savefig(OUT/'common_pair_profiles.pdf');fig.savefig(OUT/'common_pair_profiles.png',dpi=190);plt.close(fig)
lines=['# 两个新档案吸收系统：共同 Fe II 线对的不确定性','',
'这次实际计算了相同观测量 D = v(Fe II 2382) − v(Fe II 2374) 的 nuisance-profile 轮廓。两组数据都能得到数值封闭的条件轮廓；但 J004 的计算发现了此前未找到的更低气体/仪器宽度参数分支，因此不能把最初的局部标准误差直接当作可靠测量精度。','',
'结果仍然只描述指定气体结构、参数边界和误差模型下的相对线位移。仪器残余标定、像素协方差、气体模型选择与原子响应的不确定性没有被完整边缘化，本轮没有证明精细结构常数或宇宙时间律变化。','',
'## Material Passport','',
'- 模式：实际数据上的追加数值检验；原始档案模型和历史结果保持不变。',
'- 数据：J232128−105122，z=1.629，3 个气体分量/5 条线/130 像素；J004131−493611，z=2.248，6 个气体分量/3 条线/138 像素。',
'- 权重：SQUAD FITS 主数组 row2 expected fluctuation，仍为对角权重，不是完整协方差。',
'- 固定：此前选定的谱线、窗口、原始像素、掩膜、分量数和所有边界；像素积分 OS33。',
'- 独立验证：见 [独立复核](common_pair_independent_review.md) 与 `results/common_pair/review/validation.json`。',
'', '## 共同观测量的结果','',
'下表的 profile 区间是实际非线性拟合在 Δχ²=1 处的插值交点，保留了非对称性。它是条件轮廓，并未校准为具有指定覆盖率的置信区间。局部 σ 则允许 nuisance 参数在切空间中双向变化；有硬边界时，它与实际受限 profile 不必一致。','',
'| 目标 | z | D，m/s | 局部切线 σ，m/s | profile Δχ²=1 交点，m/s | profile Δχ²=3.84 交点，m/s |',
'|---|---:|---:|---:|---|---|']
for r in rows:lines.append(f'| {r["target"]} | {r["z_abs"]:.3f} | {r["D_m_s"]:.3f} | {r["conditional_sigma_m_s"]:.3f} | {fpair(interval(r,1.))} | {fpair(interval(r,3.84))} |')
lines+=['','![共同线对条件 profile](../results/common_pair/common_pair_profiles.png)','',
'图中蓝线为逐点固定 D、重新优化其他参数的实际结果；橙虚线为局部切线近似。阴影表示 Δχ²=1 的插值范围。J004 的受限非线性轮廓与无约束局部近似明显不同，参数硬边界可能贡献了这种差异；尚未将其与非线性及残差曲率的作用分别隔离。较窄的蓝色轮廓不能解释为模型独立的精度提升。','',
'J232 的最佳点有 `zero_2260` 和 `zero_2374` 两个边界活跃；J004 还涉及 `logb_1` 和 `logfwhm_1608`，以及 `zero_2374`。这些结果需要随模型边界一起使用。轮廓只是所检查分支附近的相连范围，不保证找遍所有不相连的似然区域。','',
'## J004 的新参数分支与公平重算','',
'最初从保存的 H0/H1 重新开始自由拟合，仍停在此前 H1 χ²≈104.534 的分支。固定 D≈+189 m/s 的 profile 探索进入更低分支，解除该约束后得到 χ²≈87.255、D≈−120.110 m/s。原始像素、误差和分量数都没有变化；改善来自数值搜索找到新的 nuisance 参数组合。','',
'为公平比较，随后把这组气体参数用于所有相对位移均为零的 H0 重启，并从该 H0 再启动自由 H1。下表是新 companion 中保存的有限终点比较，历史 row2 文件仍保留原值。','',
'| 目标 | 历史 H0 χ² | 历史 H1 χ² | 新 H0 χ² | 新 H1 χ² | 新 Δχ² | 新增位移参数 |',
'|---|---:|---:|---:|---:|---:|---:|']
for r in rows:
 c=r['full_comparison'];lines.append(f'| {r["target"]} | {r["original_h0_chi2"]:.6f} | {r["original_h1_chi2"]:.6f} | {c["chi2_null"]:.6f} | {c["chi2_alternative"]:.6f} | {c["delta_chi2"]:.6f} | {c["added_parameters"]} |')
lines+=['',
'J004 的最佳 D 从历史分支的约 +26 m/s 改为约 −120 m/s，本身表明仅依赖一个成功终止端点并不足以评估精度。新 H0 同样得到明显改善，因此不能把 H1 的全部改善归因于某种常数变化。','',
'此外，D=0 的 profile 仍允许其他谱线位移浮动，它不同于表中全部位移均为零的 H0。对 J004，D=0 profile 的 χ²≈87.844855，目标增加约 0.590；全部位移为零的 H0 则为约 95.228461。','',
'## 方法、数值范围与限制','',
'局部尺度通过将目标位移的加权 Jacobian 列投影到所有其他 H1 参数张成空间的正交补得到，σ_D=1000/‖j_perp‖；1000 将 km/s 换算为 m/s。列归一化 SVD 使用 10⁻⁸、10⁻¹⁰、10⁻¹² 三档阈值，没有用 reduced χ² 缩放标准误差。','',
'实际 profile 每次只固定 `shift_2382`。气体柱密度、速度、展宽、各线连续谱、零点、Gaussian LSF 宽度及其余相对位移都按原 H1 边界重新优化。主网格覆盖约 ±3 个局部尺度并受原 ±1 km/s 位移边界限制，每点至少两个起点；阈值交点再细分。J004 的改进网格另按 D 降序，从当前最低目标端点重新起步复查，不把这种复查描述为完整的连续追踪所有分支。','',
'每次优化的预算为 300 次评估；少数自由 H0/H1 重启为 400 次。所有选定 profile 端点均成功终止，两个约定阈值在两侧都已夹住。原始数值尝试全部保留；成功主要对应 ftol 停止条件，不认证全局最小值。','',
'交点用实际网格夹逼后线性插值。JSON 保留 bracket 范围，可见网格分辨率；表格多位数字便于重算，不代表天文物理精度。最优固定坐标端点与自由 H1 的目标差仅约 10⁻⁹–10⁻⁸，低于预设 10⁻⁴ 的基线一致性容差。','',
'这些数据让不同红移目标之间比较同一个谱线对成为可计算的问题；但目标环境、仪器和原子响应仍需要独立控制。这里只提供条件观测描述，没有拟合 1/ln²(t)，也没有把两个不同目标解释为已测出的宇宙时间变化。','',
'## 可复查文件','',
'- [最终拟合和 profile 汇总](../results/common_pair/refined_summary.json)。',
'- [跨目标统一字段](../results/common_pair/harmonized_measurements.json)：显式标注 `precision_measurement=false` 和 `physical_time_inference=false`。',
'- [可导出 PDF 图](../results/common_pair/common_pair_profiles.pdf)。',
'- [独立验证 JSON](../results/common_pair/review/validation.json)。',
'- 原始 profile 阶段的最后一个汇总语句出现重复 metadata-key 异常，发生在全部数值拟合写盘之后。原脚本与数值产品按原哈希保留；`common_pair_refine.py` 从这些保存结果完成新的参数分支检查和最终汇总。本报告及图由 `common_pair_report.py` 导出，不依赖那条出错的第一阶段汇总语句。','',
'```bash',
'# 从归档的第一阶段端点重现 companion 检查和导出（不覆盖历史档案结果）',
'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_clock/experiments/highz_feasibility_2026-09-20/code/common_pair_refine.py',
'python riemann_clock/experiments/highz_feasibility_2026-09-20/code/common_pair_validate.py',
'python riemann_clock/experiments/highz_feasibility_2026-09-20/code/common_pair_report.py',
'```','']
(ROOT/'reports/common_pair_results_cn.md').write_text('\n'.join(lines))
provenance=dict(input_summary_sha256=sha(OUT/'refined_summary.json'),export_script_sha256=sha(__file__),harmonized_sha256=sha(OUT/'harmonized_measurements.json'),report_sha256=sha(ROOT/'reports/common_pair_results_cn.md'),figure_pdf_sha256=sha(OUT/'common_pair_profiles.pdf'))
writepath=OUT/'export_provenance.json';writepath.write_text(json.dumps(provenance,indent=2)+'\n')
print(json.dumps(provenance,indent=2))
