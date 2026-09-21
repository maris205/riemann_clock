#!/usr/bin/env python3
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[2]
v=json.loads((ROOT/'results/followup_review/validation.json').read_text())
comparisons=v['comparisons'];all_done=v['all_comparisons_completed']
lines=['# 饱和核掩膜与经验协方差后续检验：独立数值复核','',
'## Material Passport','',
'- Task: 对固定饱和核掩膜与连续谱估计协方差的非线性检验进行独立代码与数值审查。',
'- Mode: validate；复核者没有实现被审查的拟合脚本，也没有重新运行昂贵优化。',
'- Source: 本项目归档光谱、此前 H0 拟合、连续谱控制结果，以及本轮新增输出。',
f'- Verification status: {v["verification"]}；{v["n_pass"]}/{v["n_checks"]} 项数值检查通过；' + ('三组最终比较均已复核。' if all_done else '部分最终拟合仍在运行，当前报告为阶段性结果。'),
'- Reproducibility scope: 已保存数值结果的重算为 VERIFIED；全局优化与物理解释仅为 ANALYZED。',
f'- Audit timestamp UTC: {v["timestamp_utc"]}',
'', '## 结论与范围', '',
'在本报告覆盖的输出中，未发现固定掩膜、残差/Jacobian 白化、缓存或目标函数存储存在实现错误。独立验证脚本为 `results/followup_review/verify_followup.py`；完整检查、数值误差与 SHA256 见 `results/followup_review/validation.json`。',
'',
'这些检查支持将本轮结果作为有明确条件的敏感性分析。它们不证明找到了全局最小值，不校准“发现显著性”，也不能把单个吸收系统中的相对线位移解释为精细结构常数变化或宇宙时间律。',
'', '## 已完成的配对结果', '',
'下表中的差值是两套约束非线性优化目标之差；它不同于此前固定基线、投影 nuisance 方向所得的局部分数。每组 H0/H1 使用完全相同的像素。', '',
'| 比较 | H0 | H1 | 像素数 | 额外参数 | 非线性 Δχ² | 状态 |',
'|---|---|---|---:|---:|---:|---|']
for c in comparisons:
 name=Path(c['file']).parent.name+'/'+Path(c['file']).name
 if c['status']=='PENDING':lines.append(f'| {name} | — | — | — | — | — | 运行中 |')
 else:lines.append(f'| {name} | `{c["null"]}` | `{c["alternative"]}` | {c["ndata"]} | {c["added_parameters"]} | {c["delta_chi2"]:.8f} | 已重算 |')
lines+=['','## 固定饱和核掩膜','',
'用独立的像素距离构造重建了两个掩膜，而非重复调用拟合代码中的 `binary_dilation`。同时由原物理模型重算冻结的 H0 光谱，确认阈值来源确为保存的模型。H0/H1 掩膜一致，参数变化不会移动掩膜，弱线像素没有变化，原来无效的像素没有被重新纳入。', '',
'| 掩膜 | F < 阈值 | 要求的扩张 | 实际扩张 | 删除 Fe II 2382 像素 | 删除 Fe II 2600 像素 | 保留总像素 |',
'|---|---:|---:|---:|---:|---:|---:|']
for m in v['metrics']['mask_cases']:
 counts={x['line']:x for x in m['per_line']}
 lines.append(f'| {m["case"]} | {m["threshold"]} | ±{m["requested_padding_km_s"]:g} km/s | ±{m["dilation_pixels"]} 像素 = ±{m["actual_padding_km_s"]:.6f} km/s | {counts[2382]["removed"]} | {counts[2600]["removed"]} | {m["ndata"]} |')
lines+=['',
'“±1 km/s”要求向上取整为三个原始像素；正文应明确实际约为 ±1.2 km/s。掩膜虽然在当前配对拟合中固定，但其参考光谱来自对同一批数据的先前拟合，因此它是数据辅助的探索性控制，并非盲法确认试验。删除更宽区域后偏好的降低同时可能来自信息损失，不能独自判定饱和效应就是原因。',
'', '## 经验协方差与 GLS', '',
'独立构造 C_ij = k(|index_i-index_j|)，以原光谱像素索引距离建立各谱线的协方差。用直接稠密线性方程求解验证 rᵀC⁻¹r 及 JᵀC⁻¹r；另验证 Cholesky 重建、参数缓存往返与原地修改、零相对位移时 H1 包含 H0、以及有限差分 Jacobian。前向光谱与原物理模型完全一致。', '',
f'- 连续谱直接测得 lag-one 相关：{v["metrics"]["covariance"]["raw_rho1"]:.9f}。',
f'- Bartlett taper 后用于 GLS 的 lag-one 核：{v["metrics"]["covariance"]["tapered_rho1"]:.9f}；不是未经修改的原始相关系数。',
f'- 六个协方差块均正定，最小特征值范围 {min(x["min_eigenvalue"] for x in v["metrics"]["covariance"]["blocks"]):.6f}–{max(x["min_eigenvalue"] for x in v["metrics"]["covariance"]["blocks"]):.6f}。',
f'- 检查 13 个 Jacobian 列，覆盖气体参数、连续谱及全部五个相对位移；最大相对 L2 有限差分误差 {max(x["relative_l2_error"] for x in v["metrics"]["finite_differences"]):.3g}。',
'- 三个掩膜缺口的相邻保留像素索引距离为 23、64、18，都超过十阶协方差核的支持范围；此处按原索引距离处理与在这些缺口重置的结果一致。',
f'- 主检验沿用给定误差幅度。把连续谱方差 {v["metrics"]["covariance"]["continuum_variance"]:.9f} 同样迁移到吸收区，会把所有目标函数与差值乘以其倒数；共同标量不改变最佳参数。',
'',
'协方差在 H0/H1 配对内固定，因此高斯似然中的行列式项抵消。但连续谱协方差形状是否适用于吸收区、其估计误差、跨谱段相关性和仪器校准仍未完全计入；目标差值不能直接转写为物理发现的概率。',
'', '## 优化与出处', '',
f'- 原核心 `espresso_conventional_null.py` SHA256 仍为 `{v["source_hashes"]["code/espresso_conventional_null.py"]}`。',
f'- 当前重算 {len(v["fit_replays"])} 个已完成拟合的目标值、参数、残差、Jacobian、模型、掩膜与初始参数映射。',
'- 每个结果中的冻结掩膜或协方差来源哈希均与当前归档文件一致。',
'- 对 `espresso_saturation_controls.py`、`espresso_empirical_gls.py` 和 `polish_espresso_followup.py` 均进行了源代码审查并记录哈希。继续拟合脚本保留原来达到评估次数上限的输出，重用记录的继续拟合产品，并在选择时阻止遗漏目标值更低的未终止端点；该逻辑没有发现阻断性缺陷。',
'- 已完成的最终配对均检查成功状态、五个额外参数及最佳已终止初值选择；另检查是否存在目标值更低但尚未成功终止的已保存尝试。',
'- `optimizer_success=True` 在本轮主要对应 `ftol` 停止条件，而非已达到很小投影梯度。多个参数仍在边界附近，气体成分方向存在退化。准确表述为“在所测试初值中、按指定容差终止的最佳拟合”，不应声称已证明全局收敛。',
'', '## 统计推断边界扫描（11/11）', '',
'| 项目 | 本任务中的核查结论 |',
'|---|---|',
'| Simpson 悖论 | 不把不同掩膜、不同信息量的目标差值当作同一统计量直接合并；按谱线保存输出。 |',
'| 生态谬误 | 单一吸收系统的谱线不一致不能推广为宇宙时间的演化关系。 |',
'| 选择偏差 | 已采用发表的系统/窗口；新掩膜由此前拟合辅助定义，须保持探索性标记。 |',
'| Collider 偏差 | 无因果效应估计；同数据辅助定义的条件筛选不宣称独立。 |',
'| 基准率忽略 | 未估计新物理的后验概率；不把目标差值当成该概率。 |',
'| 均值回归 | 删除问题区域后的统计量降低不构成机制证据；独立仪器确认需另检验。 |',
'| 幸存者偏差 | 保留原始和未成功尝试；最终选择审计包括是否有更低的未终止目标。 |',
'| 多重比较/别处效应 | 两个阈值和多个协方差敏感性模型属于探索；没有给出已校准的发现显著性。 |',
'| 分叉路径 | 所有设置、起点和输出保持可追踪；不把事后分析写成预注册检验。 |',
'| 相关与因果 | 相对线位移不能单独归因于 α 演化、宇宙时钟或某个仪器机制。 |',
'| 反向因果 | 单红移截面没有辨识时间方向，更没有检验 1/ln²(t)。 |',
'', '## 重现', '',
'```bash',
'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_clock/experiments/highz_feasibility_2026-09-20/results/followup_review/verify_followup.py',
'python riemann_clock/experiments/highz_feasibility_2026-09-20/results/followup_review/write_review.py',
'```','']
report=ROOT/'reports/espresso_followup_independent_review.md';report.write_text('\n'.join(lines))
print(report)
