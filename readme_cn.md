# Riemann Clock：有限资源下的黎曼零点估计

[English](readme.md)

本项目重新整理“物理测量的黎曼零点精度有限，其可达高度可能随宇宙时间增加”的设想。**数学零点本身保持不变；拟研究的是指定物理系统如何估计它，以及估计能力是否存在额外的宇宙时间依赖。**

新的论文、数据、代码、图和报告全部在 `riemann_clock`。原 `Cosmic-Chaos-Alpha` 项目保持原样。

## 先看这些

- [原假说论文：35 页英文单栏 PDF](paper/riemann_clock.pdf)
- [新增光谱分析专稿 PDF](paper_spectroscopy/feii_relative_frequencies.pdf)
- [中文核查与修订报告](reports/revision_report_cn.md)
- [独立数值复核](reports/numerical_validation.md)
- [实验来源核查](reports/experimental_sources_audit.md)
- [引文核查](reports/citation_checks.md)
- [λ-cosmos 综合稿的补充对照与修订建议](reports/lambda_cosmos_comparison_cn.md)（单独核查，尚未合入正文）
- [高红移光谱：真实数据先导与十年实验方案](experiments/highz_feasibility_2026-09-20/experiment_protocol_cn.md)（完整方案作为支持材料；实测分析与展望已合入正文）
- [已有真实数据分析与未来观测展望](experiments/highz_feasibility_2026-09-20/reports/existing_data_results_and_outlook_cn.md)（实际通量测量和配准诊断已加入论文）

论文是待作者审阅的工作稿。没有声称已经完成外部同行评审，也没有代作者确认资助、利益冲突等投稿声明。

## 新增：独立的光谱分析专稿

[新单栏论文](paper_spectroscopy/feii_relative_frequencies.pdf) *Instrument and Gas-Model Sensitivity of Fe II Relative-Frequency Tests in Archival Quasar Spectra* 将天文学部分整理为独立研究，集中讨论气体结构、协方差、原生曝光、相邻级次、共同线对和跨年龄可识别性。[LaTeX 源码](paper_spectroscopy/main.tex) 与原来的 **35 页假说稿**分开保存；原稿的实验室分析、理论来由和资源热力图保留。新专稿第5.4节简述预测实验，附录A给出时间函数、独立验证、协方差、失败判据、条件功效和仪器要求；明确这仍是尚未取得新观测、未冻结物理幅度的设计。[完整设计记录](experiments/prediction_test_2026-09-21/prediction_protocol_cn.md)一并保存。

正值、零质心的非对称响应控制下，Fe II 2600 级次差仍约 **−59 m/s**，但尚未定位其仪器、提取或模板原因。共同线对的条件轮廓还找到了更低的气体参数分支，说明此前“成功终止”的数值解并非最终精度保证。下面旧阶段的数值按历史记录保留，最新解释以新专稿和本轮报告为准。跨年龄曲线比较属于指定振幅下的可识别性计算，没有测出宇宙演化振幅或常数变化。

- [本轮完整结果与论文整理说明](experiments/highz_feasibility_2026-09-20/reports/publication_followup_results_cn.md)
- [级次响应模型](experiments/highz_feasibility_2026-09-20/reports/order_response_results_cn.md)与[残差、迹线控制](experiments/highz_feasibility_2026-09-20/reports/order_response_residuals_cn.md)
- [共同 2382/2374 线对及条件轮廓](experiments/highz_feasibility_2026-09-20/reports/common_pair_results_cn.md)
- [跨年龄可识别性与函数区分能力](experiments/highz_feasibility_2026-09-20/reports/cross_age_identifiability_cn.md)
- [历史标定资料和后续需求](experiments/highz_feasibility_2026-09-20/reports/calibration_readiness_cn.md)

在本项目目录用已有结果与图片编译新稿：

```bash
python code/build_spectroscopy_paper.py
```

该命令不改写原假说论文。两份论文均未声称获得外部同行评审或期刊接收。

## 前一轮：保留在原假说稿中的已有数据扩展

已分析**全部17次ESPRESSO原生曝光**，新增读取**32份UVES光谱、36个吸收系统、324个谱线窗口**。6个系统通过自动多谱线质量筛查，其中两个完成共享气体的物理拟合。这些是历史观测的新分析，不是2026年新取得的光子，也不是6个精密变化测量。

由连续谱相关构造的协方差模型与强线核心剔除联合使用后，五个位移的改善为**Δχ²=23.69、6.70**。宽遮罩同时移除了指定两线共同偏移方向约84%的局部信息，因此改善减弱不能单独定位饱和原因。固定40分量气体替代给出18.87，并显示明显初值依赖。

一个具体的新线索是：**同一Fe II 2600在相邻光栅级次中相差约60 m/s，放开Gaussian仪器宽度后仍约58 m/s**。这值得优先检查仪器、提取、标定和模板近似，但还未定位原有全部偏移的原因。两个新目标的常规模型已能较好拟合；使用档案推荐误差数组，额外位移只改善6.16（4参数）和5.50（2参数）。严格UVES续算也如实保留达到次数上限、未满足一阶驻点要求的状态。

**目前是把“是否超出常规解释”这一关做得更扎实，尚未发现物理常数变化，也未拟合出对数平方反比时间规律。**

- [本轮完整结果与解释](experiments/highz_feasibility_2026-09-20/reports/available_data_campaign_results_cn.md)
- [联合控制及信息损失图](figures/joint_controls_information.png)
- [17次曝光、批次和级次诊断图](experiments/highz_feasibility_2026-09-20/results/exposures/exposure_consistency.png)
- [36个系统的谱线质量热力图](experiments/highz_feasibility_2026-09-20/results/archive_expansion/archive_line_readiness.png)
- [两个新目标的实际物理拟合](experiments/highz_feasibility_2026-09-20/reports/archive_profile_results_cn.md)
- [最终PDF与交付一致性检查](reports/available_data_delivery_verification.md)

本轮方法、比较表和联合信息图已合入单栏论文。下文保留此前各阶段说明；原资源热力图及实验室图表保留。

## 保留的核心假设

设指定测量协议中存在一个额外的标准差分量：

```math
\delta_{\mathrm c}(t)=\delta_{\mathrm c0}
\bigg[\frac{\ln(t_0/t_{\ast})}{\ln(t/t_{\ast})}\bigg]^2.
```

这里对数的自变量无量纲，振幅可以为零。形式来自已发表的[非自治二次映射与黎曼零点数值对应论文](https://doi.org/10.3390/mca31050193)的启发。但是“映射迭代次数 → 宇宙年龄 → 可测噪声或精度”的两步连接仍需物理机制，不能当作原论文已证明的结论。

取 138 亿年和普朗克时间作为示例尺度，假设分量的现今相对变化率约为 **−1.03 × 10⁻¹²/年**，半年约 **−5.17 × 10⁻¹³**。这是条件计算，没有用实验拟合出振幅或发现宇宙老化。

## 本次数据更新

2021 年离子阱原始补充表转录出 **269 条估计**，涉及 80 个零点、四个驱动设置。新版按独立计算的数学参考值重算残差，正确解释括号误差，并保留原表中的参考值差异。部分误差条增大，但没有证实共同的“第 80 个零点物理极限”。

找到一篇 **2026 年 7 月 1 日正式发表**的核自旋研究，加入公开的 **110 个时间扫描点、18 个逆温度扫描点**和论文报告的五个根估计。预印本在 2025 年 11 月已出现，不能称为最近半年首次报告。硬件实验研究前五个零点附近；第万亿个零点附近属于数值模拟。

这些不同平台的数据不能拼成“宇宙时间变长，能看到的零点就更多”的证据序列。最近实验让协议研究更丰富；宇宙时间部分仍是假设。

## 真实天文光谱分析

新版还加入两份 SQUAD 真实合并光谱的分析，吸收红移分别为 1.1508 和 3.025。按 NIST Ritz 波长计算有限窗口的吸收积分、质心及稳健性，并用真实 Fe II 2374 轮廓配准 2586、2600，检查模型残差。新图和数据表已加入单栏论文；原有四张图保留。

HE 的吸收质心条件统计误差为 18–117 m/s，但人为设置的 1% 连续谱扰动可改变质心 99–439 m/s，边界变化也会改变包含的气体分量。高红移样例缺少已确认干净的多条差分跃迁，因此没有报告宇宙演化系数。未来设计据此优先处理线组、混叠、饱和和校准。

- [真实数据诊断图](figures/highz_observed_profiles.png)
- [真实通量的经验模板拟合](figures/highz_observed_registration.png)
- [实测数据与代码目录](experiments/highz_feasibility_2026-09-20/readme.md)

此前的人工谱线注入与 24／60 视线功效预测继续单独标为模拟；它们不是这两份实际光谱已经达到的宇宙参数精度。

## 重现

在本目录运行，依赖见 `requirements.txt`；还需 `pdftotext`，编译论文另需 `pdflatex`、`bibtex`、REVTeX 4.2。

```bash
python code/extract_public_experiments.py
python code/analyze_zero_precision.py
python code/plot_recent_experiment.py
python code/plot_resource_heatmap.py
python code/validate_analysis.py
python code/build_paper.py
```

重现新增真实光谱分析（依赖该实验目录中的 `requirements.txt`，使用已有缓存无需联网）：

```bash
python experiments/highz_feasibility_2026-09-20/code/measure_observed_profiles.py
python experiments/highz_feasibility_2026-09-20/code/register_observed_profiles.py
python experiments/highz_feasibility_2026-09-20/code/design_from_observed_data.py
python experiments/highz_feasibility_2026-09-20/code/prepare_highz_manuscript.py
python code/build_paper.py
```

已有下载文件可离线重现。`data/raw` 保留公开来源和作者提供的处理后数据，并非完整原始采样记录；`data/processed` 是转录与诊断；`results` 是数值摘要和明确标注的条件情景；`figures` 包含四张原有新版图和两张新增天文分析图。天文来源与完整结果保存在上述实验目录。旧项目中人为构造“误差悬崖”的模拟图没有作为实验资料使用。

新增的[资源—精度热力图](figures/resource_precision_heatmap.png)保留了渐变配色和星形标记。横轴是重复测量次数，纵轴是假设的精度底，颜色表示模型中零点高度包络的十进制对数；星标 A–C 是示例参数，虚线区分统计项与精度底主导的区域。

## 前一阶段：常规吸收物理基线检验

已补做真实 ESPRESSO 光谱的联合 Voigt 检验：六条 Fe II 谱线、2931个有效像素、45个共享气体分量。增加五个相对位移后，基准比较改善 **Δχ²=51.71**，主要涉及2382和2600两条强线。数项常规模型/误差对照仍保留条件性偏好；去掉两条强线后改善降至7.09，对应三个额外参数。

准确结论是：**指定模型下存在值得复核的相对谱线不一致，尚未建立超出常规解释的物理变化。** 这不是精细结构常数测量，也没有据此拟合宇宙时间的对数平方反比规律。

- [本轮完整核查与结果](experiments/highz_feasibility_2026-09-20/reports/conventional_null_results_cn.md)
- [实测谱线与常规模型图](figures/espresso_conventional_null.png)
- [控制比较图](figures/espresso_null_robustness.png)
- [代码、原始资料与复现说明](experiments/highz_feasibility_2026-09-20/readme.md#espresso-conventional-null-test)

结果已单独成节加入单栏论文。既有热力图与实验室分析保留。

## 前一阶段：分别处理噪声、强线核心和独立仪器

在 25 个不重叠连续谱区间的 13,750 个像素中，实测邻像素相关约 **0.373**。以实测相关形状构造固定协方差后重新进行非线性拟合，五个相对位移的改善从 51.71 降到 **29.24**。分别遮掉两条强线最深核心、或更宽核心并外扩三像素，改善为 **41.66、13.26**。这些是分别进行的条件比较，不能当成已同时排除全部常规因素。

UVES 同一吸收系统的共享气体拟合也已开展，但其仪器线宽、触边参数与局部解限制了独立验证。当前仍是“相对谱线不一致值得追踪”，没有建立物理常数变化，也没有拟合出宇宙时间规律。

- [本轮结果与解释](experiments/highz_feasibility_2026-09-20/reports/conventional_followup_results_cn.md)
- [连续谱相关及非线性重拟合图](figures/espresso_noise_controls.png)
- [UVES 交叉检验及限制](experiments/highz_feasibility_2026-09-20/reports/cross_instrument_results_cn.md)

2026 年 7 月的 [HARPERFECT](https://arxiv.org/abs/2607.06809) 提供相关的新提取方法，处理的仍是历史 HARPS 曝光；本次未取得科学数组与配套分辨率矩阵，因此没有声称已分析这些产品。
