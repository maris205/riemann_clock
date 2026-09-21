# UVES 独立仪器复核：条件性同方向线索，尚未排除常规解释

本轮实际拟合了独立的 UVES 科学光谱。所用三条避开已知混线的 Fe II 谱线，在这套物理模型下给出与 ESPRESSO 相同方向的相对位移；但气体模型存在不同局部解，本轮常规模型端点均未达到预设收敛准则，且参考线的拟合仪器线宽触及预设边界。因此它是值得保存的**条件性相容线索**，不能称为已经确立的独立测量，更不能据此宣布物理常数变化。

## 数据和模型

使用 [UVES SQUAD DR1](https://github.com/MTMurphy77/UVES_SQUAD_DR1) 的 HE0515−4414 合并光谱：51 次曝光，181652.0812 秒；输入记录日期为 1999-12-14 至 2009-02-11。这不是 2026 年新观测。真空日心波长网格约 1.3 km/s；直接使用原始合并像素，没有再次重采样观测通量。

选取 Fe II 2374、2382、2600，共 511 个有效像素（156、178、177）。2586 有已知额外吸收，本轮主比较没有使用；2344 也排除。适配代码保留了 2586 遮罩接口，但不得把未执行的四线拟合写成已完成。

常规模型包含同位素分裂、自然展宽、45 个共享气体分量的柱密度/速度/Doppler 宽度、像素积分、每条线的连续谱归一化和斜率、零点修正、对称高斯仪器线形。45 分量架构由 ESPRESSO 已发表模型提供；所有这些气体数值在 UVES 中重拟合。原子波长、同位素结构和振子强度也使用与 ESPRESSO 相同的数据表。两项分析共享这些物理输入与模型架构；不同仪器提供独立的观测光子数据。

H0 固定实验室跃迁关系；H1 在完全相同的干扰参数模型上，仅增加 2382 和 2600 相对 2374 的两个速度偏移。常规拟合有 147 个参数，偏移模型有 149 个；不透明度敏感性分支分别增加两个参数。高斯 FWHM 的允许范围为档案名义值的 0.65–1.25 倍。±10% 是明确指定的敏感性范围，不是从实验室测量推导出的先验误差分布。

## 保存结果

以下 χ² 差值是**有限优化和指定模型下的目标函数差**。三组比较均有尚未满足收敛准则的端点，因此差值均属暂定；没有换算为“几 σ”，没有使用名义自由度计算发现概率。

| 分支 | H0 χ² | H1 χ² | 差值 | 2382−2374 (m/s) | 2600−2374 (m/s) |
|---|---:|---:|---:|---:|---:|
| 统计误差 | 417.841 | 406.847 | 10.993 | -139.2 | -131.7 |
| 档案建议的 expected-fluctuation 误差 | 398.118 | 388.135 | 9.982 | -140.2 | -132.3 |
| 同上，并允许强线不透明度各变动 ±10% | 396.724 | 386.570 | 10.154 | -139.8 | -130.9 |

每个比较选择同一分支内保存的最低 χ² 端点；候选文件全部保留在 `results/cross_instrument/summary.json`。优化终止情况如下。函数评估数是该次续算的次数，不是累计总次数；续算源文件记录在各 JSON 的 `input_start_file` 中。

| 端点文件 | 终止情况 | 本次函数评估次数 nfev | 最优性指标 | 数值秩/参数数 | 活动边界数 |
|---|---|---:|---:|---:|---:|
| `uves_refined_null_continued` | 达到迭代上限 | 500 | 0.1485 | 147/147 | 11 |
| `uves_refined_alternative` | 达到终止容差 | 554 | 0.001895 | 149/149 | 12 |
| `uves_expected_null_continued` | 达到迭代上限 | 500 | 0.387 | 147/147 | 12 |
| `uves_expected_alternative_continued` | 达到终止容差 | 385 | 0.04115 | 149/149 | 12 |
| `uves_opacity_null` | 达到迭代上限 | 350 | 1.322 | 149/149 | 12 |
| `uves_opacity_alternative` | 达到迭代上限 | 350 | 0.1589 | 151/151 | 12 |

初始统计误差 H0 分支停在 χ² = 427.119；从 H1 的气体参数重新启动 H0 后得到更低的解。这个变化本身说明不能只挑一个初始化状态就宣称显著性。相对位移的优化解与仪器、气体参数的局部解存在耦合。

## 关键的常规解释与限制

档案的名义分辨率 R≈53696、弧灯分辨率 R≈70723，都不是独立测得的类星体有效线形。直接固定名义高斯线宽时 H0 χ²≈1399；允许拟合仪器线宽后目标函数降低数百，证明仪器线形假设会显著影响此项检验。

最终统计误差 H1 的高斯 FWHM 分别为 3.6290、3.6969、3.7390 km/s。2374 的宽度比例达到预设下限 0.65（相对于名义宽度），对应有效 R≈82609。不能把这个边界上的拟合值当作仪器标定测量。气体速度和宽度也有活动边界，因此局部协方差误差不能涵盖全部不确定性。

该合并光谱并非简单“未经修正”。原始 UPL 文件的 138 条谱臂记录中有 95 条非零 VSHT，保存了逐曝光波长位移/斜率修正。官方 [UVES_popler 头文件](https://raw.githubusercontent.com/MTMurphy77/UVES_popler/master/UVES_popler.h) 和 [波长计算代码](https://raw.githubusercontent.com/MTMurphy77/UVES_popler/master/UVES_wpol.c) 说明其含义。不过，本次没有逐曝光重做标定，也没有确认它与 Kotuš 等精密分析的最终产品完全等价；剩余标定误差和完整像素协方差尚未重建。

档案建议用 expected-fluctuation 数组做 χ² 分析，已实际执行这一分支。它相对于统计误差的中位比例在三条线分别约为 1.0823、1.0071、1.0034，尤其改变参考线权重，但仍未提供完整像素协方差。

两个仪器观测同一吸收系统。即便一个相对位移模式跨仪器出现，也可能来自共享的气体结构、原子强度、同位素组成或模型假设；这不是两个独立宇宙年代的数据点。没有拟合 α、黎曼零点截断参数或 1/ln²(t)。

## 数值与来源核查

独立复核脚本为 `code/cross_instrument_validate.py`。在本报告生成时，机器报告为 PASS，310/310 项检查通过；最终审计覆盖的具体端点和 SHA-256 均在 `results/cross_instrument/validation.json`，检查通过只表示实现/保存数组可复算，不表示物理推断已获验证。

早期每像素 7 子点的积分对 χ² 有约 1 的影响，所以仅保留为探索记录；本表统一采用 21 子点，并在最终参数上与更密网格独立比较。初始分支/中间检查点不能替代本表指定的最终端点。

[HARPS 2021 原文](https://arxiv.org/html/2008.10619v2) 明确声明提供 LFC/ThAr 合并光谱和模型，但本次未能从当前出版商入口取得可验证的附件数组；不得表述为作者未公开数据。[2026 年 7 月 HARPERFECT](https://arxiv.org/html/2607.06809v1) 是同一历史 52.5 小时 HARPS 观测的新提取方法，有望改善线形和协方差处理，但此次检索未找到可下载的科学数组/分辨率矩阵。详见 `reports/cross_instrument_harps_availability.md`。

## 复现和产物

代码：`code/cross_instrument_uves.py`、`cross_instrument_validate.py`、`cross_instrument_summarize.py`、`cross_instrument_report.py`。

所有模型参数、所用像素、噪声数组、模型通量、Jacobian、局部协方差和终止状态保留在 `results/cross_instrument/`。最终参数数值复核不需重新优化：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_clock/experiments/highz_feasibility_2026-09-20/code/cross_instrument_validate.py
python riemann_clock/experiments/highz_feasibility_2026-09-20/code/cross_instrument_report.py
```

示例拟合命令（具体续算起点见对应 JSON）：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_clock/experiments/highz_feasibility_2026-09-20/code/cross_instrument_uves.py --name uves_example --free --flex-lsf --error-row 2 --oversample 21 --max-nfev 500 --start riemann_clock/experiments/highz_feasibility_2026-09-20/results/cross_instrument/uves_expected_alternative.json
```

图：`results/cross_instrument/uves_profiles.pdf` 和 `.png`；图中为收敛的统计误差 H1 对实际 UVES 通量的描述，不表示已经发现常规理论以外的变化。
