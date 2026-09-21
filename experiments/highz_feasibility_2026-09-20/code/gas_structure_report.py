#!/usr/bin/env python3
"""Generate a scoped report from completed architecture-control results."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/gas_structure'


def main():
    c=json.loads((OUT/'comparison.json').read_text())
    a=json.loads((OUT/(c['null']+'.json')).read_text())
    b=json.loads((OUT/(c['alternative']+'.json')).read_text())
    base=c['baseline45']
    old=json.loads((ROOT/'results/empirical_gls'/(base['alternative']+'.json')).read_text())
    lines=[
      '# Fe II 气体分量结构敏感性结果',
      '',
      '本轮新增了一个实际重拟合的气体结构，而不是只改变原 45 分量模型的起点。本方案在已知此前 45 分量结果后设计，属于探索性控制；合并规则在查看本轮新位移前冻结。按事先固定的发表参数分组合并规则，得到 40 分量模型；两个假设在这同一个较简约架构内部比较。这里没有拟合宇宙时间规律。',
      '',
      '## 核心数值',
      '',
      '| 气体结构 | H0 参数数 | H1 参数数 | H0 目标函数 | H1 目标函数 | 增加五个位移的 Δχ² |',
      '|---|---:|---:|---:|---:|---:|',
      f'| 原 45 分量、连续谱实测 GLS 核 | 147 | 152 | {base["null_chi2"]:.6f} | {base["alternative_chi2"]:.6f} | {base["delta_chi2"]:.6f} |',
      f'| 固定合并后 40 分量、相同 GLS 核 | {c["null_parameters"]} | {c["alternative_parameters"]} | {c["null_chi2"]:.6f} | {c["alternative_chi2"]:.6f} | {c["delta_chi2"]:.6f} |',
      '',
      '所有行均使用相同的六条跃迁、2931 个有效像素、原子数据、Gaussian 仪器轮廓和固定的连续谱 ACF 衰减协方差核。连续谱实测方差因子没有用于重新缩放本表。各行内部的 H0/H1 使用相同气体及连续谱自由度，H1 增加五个以 Fe II 2374 为锚点的相对位移。跨行气体结构及局部参数边界不同；本表不把跨行目标函数差当作普通嵌套检验。',
      '',
      '| 跃迁 | 原 45 分量 GLS 位移 (m/s) | 40 分量 GLS 位移 (m/s) |',
      '|---|---:|---:|',
    ]
    for k in ['2260','2344','2382','2586','2600']:
      lines.append(f'| Fe II {k} | {old["shifts_m_s"][k]:+.3f} | {b["shifts_m_s"][k]:+.3f} |')
    lines += ['', f'这一固定简约结构下，增加位移的目标函数改善从 {base["delta_chi2"]:.2f} 降至 {c["delta_chi2"]:.2f}；两条强线位移仍为负，但其幅度相对原架构各减少约 22 m/s。位移偏好并未完全局限于原 45 分量参数化，其强度和数值也确实对气体结构敏感。40 分量的绝对目标函数更高，本检验没有证明它比原架构更真实或更优。', '', '这些是指定结构与固定误差模型下的最佳局部拟合参数，不能直接称为实验测得的物理常数漂移。单一吸收系统也不能识别 1/ln²(t/t*) 或其他宇宙时间规律。', '', '## 优化记录', '', '| 起点/终点 | χ² | 成功终止 | 函数评估次数 | optimality |', '|---|---:|:---:|---:|---:|']
    for r in c['all_final_attempts']:
      lines.append(f'| `{r["name"]}` | {r["chi2"]:.6f} | {r["optimizer_success"]} | {r["nfev"]} | {r["optimality"]:.6g} |')
    lines += ['', f'本表最终选择 H0=`{a["name"]}`、H1=`{b["name"]}`。H0 成功终止：{a["optimizer_success"]}；H1 成功终止：{b["optimizer_success"]}。求解器成功终止只表示满足其局部数值停止规则，不认证全局最小值。所有未收敛的原始尝试（若存在）仍保留在结果目录。', '', 'H1 的两个起点分别停在 2236.290295 和 2249.241251，相差约 12.95。这表明局部极小值依赖不可忽略；所选结果是本轮尝试中最好的已收敛终点，不是已认证的全局最优值。', '', f'H0 接近边界的参数：`{", ".join(a["active_bounds"]) or "无"}`。', '', f'H1 接近边界的参数：`{", ".join(b["active_bounds"]) or "无"}`。', '', '## 结构如何改变', '', '按发表的 45 分量初始速度排序，对相邻间距 <1.025 km/s 的分量连组。非单成员组为原 0 起算索引 `[13,14,15]`、`[18,19]`、`[31,32]` 和 `[35,36]`。保持总柱密度、N 加权中心和 Doppler 核二阶矩以确定初值，再重新拟合全部气体参数。成员关系及边界在优化过程中不变。', '', '每个合并后分量仍有独立的柱密度、速度和 Doppler 宽度。边界围绕由发表参数算出的固定合并初值，速度范围仍为 ±2 km/s。该架构相对原模型少了 15 个气体自由参数，但没有改变每条谱线的连续谱自由度或 H1 的位移参数数目。', '', '半个仪器 FWHM 的间隔规则只用于定义一个可复算的敏感性架构，不证明其中的真实结构不可观测。高信噪比谱线的形状可以约束部分亚分辨率结构，合并后可能损失这类信息。', '', '## 结论边界', '', '本轮能够检验“位移偏好是否完全依赖原来的 45 分量参数化”；不能穷尽气体速度分布、部分覆盖、非 Gaussian 仪器响应、不同曝光的标定误差或原子数据不确定性。所有宇宙时间或基本常数解释仍需独立的物理响应关系及跨红移样本。不得把架构控制中的 Δχ² 直接换算成新的物理发现显著性。', '', '## 可复算文件', '', '- 固定方案：`reports/gas_structure_protocol_cn.md` 与 `results/gas_structure/design.json`。', '- 拟合实现：`code/gas_structure_refit.py`；该脚本不修改原先验证过的物理前向模型。', '- 所有结果、参数、残差和解析 Jacobian：`results/gas_structure/`。', '- 汇总：`results/gas_structure/comparison.json`。', '- 独立验证：`results/gas_structure/independent_validation.json` 与 `reports/gas_structure_independent_review.md`（按其实际记录判断完成状态）。', '']
    report=ROOT/'reports/gas_structure_results_cn.md'
    report.write_text('\n'.join(lines))
    print(report)

if __name__=='__main__': main()
