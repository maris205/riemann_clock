#!/usr/bin/env python3
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[2]
v=json.loads((ROOT/'results/joint_controls_review/validation.json').read_text())
lines=['# 饱和区域剔除与经验 GLS 联合检验：独立复核', '',
'## Material Passport', '',
'- Mode: validate；独立代码审查和已保存数值结果重算。',
'- Material: 新增 `code/joint_controls_fit.py`、此前验证的物理模型/冻结掩膜、固定连续谱协方差控制及本轮配对输出。',
'- Ownership: 复核者没有实现或修改被审查的拟合脚本；未重复昂贵优化。',
f'- Verification status: {v["verification"]}；{v["n_pass"]}/{v["n_checks"]} 项检查通过。',
f'- Final comparisons available: {v["all_comparisons_completed"]}；selected optimizer termination flags: {v["all_selected_optimizers_terminated"]}。',
'- Scope: 保存结果和代数重算为 VERIFIED；全局优化、噪声迁移的物理有效性及宇宙时间解释不在数值认证范围内。',
f'- UTC: {v["timestamp_utc"]}', '',
'## 结论', '',
'所覆盖结果未发现阻断性的实现错误。联合模型先固定剔除强吸收区域，再在保留的原始像素索引上重新构造协方差并白化残差/Jacobian。H0/H1 的像素、协方差与气体/连续谱自由度一致，额外自由度仅为五个相对谱线位移。', '',
'这是对两项常规效应同时施加控制后的新一组约束非线性拟合。它不能由之前分别改变掩膜和协方差的两项 Δχ² 相加或相减得到；不同数据集或不同权重之间的差值也不是同一个似然比检验。', '',
'## 最终配对重算', '',
'| 情景 | H0 结果 | H1 结果 | 保留像素 | 增加参数 | 非线性 Δχ² | 已选优化均终止 |',
'|---|---|---|---:|---:|---:|---|']
for x in v['comparisons']:
 if x['status']=='PENDING':lines.append(f'| {x["case"]} | — | — | — | — | — | 运行中 |')
 else:lines.append(f'| {x["case"]} | `{x["null"]}` | `{x["alternative"]}` | {x["ndata"]} | {x["added_parameters"]} | {x["delta_chi2"]:.8f} | {x["selected_both_terminated"]} |')
if not v['all_comparisons_completed']:lines+=['','本报告仍为阶段性快照；最终选择待拟合完成后再审计。']
lines+=['','## 独立验证的实现细节','',
'1. 以冻结 H0 模型的阈值中心与像素距离独立重建掩膜，没有复用拟合脚本的形态扩张函数。原来无效像素没有重新纳入，弱线像素不变，改变拟合参数不会移动掩膜。',
'2. 以原始像素索引逐行构造协方差，验证它恰为原协方差在保留像素上的主子矩阵。每个 H0/H1 协方差相同且正定。',
'3. 通过独立稠密方程求解检验残差二次型及梯度，验证对白化后残差和 Jacobian 的一致处理、参数缓存的重复使用/原地改变/往返，以及前向光谱不变。',
'4. 两个掩膜分别检验 13 个有限差分 Jacobian 列，覆盖所有参数族和五个位移。',
'5. 将 H0 参数加上五个零位移嵌入 H1，验证残差和 nuisance Jacobian 一致；最终已选 H0 端点也单独进行嵌入检查。',
'6. 所有完成输出逐一重算目标函数、残差、Jacobian、模型、像素、参数边界、梯度最优性、初始参数映射和源文件哈希。',
'7. 最终选择检查覆盖全部记录端点，包括达到评估上限的端点；同一假设若存在更低目标的未终止端点，不能静默忽略。', '',
'| 掩膜 | 阈值 | 要求扩张 km/s | 实际扩张 km/s | 原始像素扩张数 | 保留像素 |',
'|---|---:|---:|---:|---:|---:|']
for x in v['cases']:lines.append(f'| {x["case"]} | {x["threshold"]} | ±{x["requested_padding_km_s"]:g} | ±{x["effective_padding_km_s"]:.6f} | ±{x["padding_native_pixels"]} | {x["ndata"]} |')
allgaps=[(x['case'],g) for x in v['cases'] for g in x['gaps']]
short=[(case,g) for case,g in allgaps if g['within_kernel']]
max_fd=max(d['relative_l2_error'] for c in v['cases'] for d in c['finite_differences'])
min_eig=min(d['min_eigenvalue'] for c in v['cases'] for d in c['covariance_blocks'])
lines+=['',f'- 有限差分最大相对 L2 误差：{max_fd:.3g}。',f'- 协方差块最小特征值：{min_eig:.6f}。',
f'- 两个掩膜合计 {len(allgaps)} 个保留像素缺口，其中 {len(short)} 个缺口仍处于十阶协方差支持范围。',
f'- 连续谱原始 lag-one 相关为 {v["raw_rho1"]:.9f}；经 Bartlett taper 后使用 {v["tapered_kernel"][1]:.9f}。',
f'- 给定误差的幅度保持原值；连续谱方差 {v["continuum_variance"]:.9f} 没有迁移到主权重中。', '']
if short:
 lines+=['新掩膜中的短缺口仍保留非零相关，验证没有把缺口两侧误当作紧邻，也没有把本应保留的跨缺口相关强制清零。', '',
'| 情景 | 谱线 | 缺口两侧原索引距离 | 协方差元素 |', '|---|---:|---:|---:|']
 for case,g in short:lines.append(f'| {case} | {g["line"]} | {g["native_distance"]} | {g["covariance"]:.9g} |')
lines+=['','## 优化与解释边界','',
'初始 H0 对两个来源只比较初始目标值，并优化其中较低者；之后再运行 H1 和由 H1 回到 H0 的交叉起点。因此不应写成“两个原始 H0 起点均完成独立优化”。', '',
'最终选择取每个假设所有记录端点中的最低目标值，并公开其终止标记。验证另确认已选 H1 目标不大于 H0，且没有遗漏更低的受限次数端点。脚本保留初始尝试与继续拟合，不覆盖原始来源。', '',
'优化成功仍仅表示满足指定的停止容差。多气体成分的退化、边界参数、协方差迁移的不确定性、仪器校准和事后选择都没有因这些数值检查而消失。', '',
'掩膜参考模型来自同一批观测的此前拟合，因此属于数据辅助的探索性控制。更宽的剔除会改变信息量，目标偏好减弱不能单独证明饱和是唯一原因；单一吸收系统也无法辨识 α 的宇宙时间变化或 1/ln²(t)。', '',
'## 统计边界扫描（11/11）','',
'| 项目 | 核查范围 |', '|---|---|',
'| Simpson 悖论 | 分开呈现两个掩膜；不把不同选择/权重的统计量混合。 |',
'| 生态谬误 | 不从一个系统推断整个宇宙的演化。 |',
'| 选择偏差 | 固定掩膜仍由同数据此前拟合辅助确定，保持探索性标记。 |',
'| Collider 偏差 | 无因果效应估计；筛选不是独立的物理确认。 |',
'| 基准率忽略 | 不把目标差值解释为新物理的后验概率。 |',
'| 均值回归 | 改变筛选后的统计量降低不独自辨识原因。 |',
'| 幸存者偏差 | 所有完成及受限次数端点保留，并审计最低目标选择。 |',
'| 多重比较 | 多个模型控制未被转写为已校准的发现显著性。 |',
'| 分叉路径 | 明示事后敏感性分析和所有已运行分支。 |',
'| 相关与因果 | 相对线位移不直接等于 α 演化或某个仪器效应。 |',
'| 反向因果 | 单红移截面没有辨识时间方向。 |', '',
'## 产物与重现','',
f'- 已重算拟合输出：{len(v["fit_replays"])} 个。',
'- `results/joint_controls_review/validation.json`：全部检查、各输出和源文件 SHA256、掩膜/协方差/有限差分度量。',
'- `code/joint_controls_validate.py`：独立验证脚本，未重复优化。',
'- 旧的 followup 审计、原掩膜比较和单独 GLS 比较在验证前后的哈希一致。', '',
'```bash',
'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_clock/experiments/highz_feasibility_2026-09-20/code/joint_controls_validate.py',
'python riemann_clock/experiments/highz_feasibility_2026-09-20/results/joint_controls_review/write_report.py',
'```','']
path=ROOT/'reports/joint_controls_independent_review.md';path.write_text('\n'.join(lines));print(path)
