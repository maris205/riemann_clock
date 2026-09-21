# 实际物理模拟实验与近半年更新核查

核查截至 2026-09-20；主要时间窗 2026-03-20 至 2026-09-20。本报告只把可核验的原始论文、作者数据和出版信息作为实验依据。没有发送邮件索取数据，没有执行下载的上游代码，没有从图片臆造点。

## 旧稿的 USTC 引文与误差需要更正

正确文献：Ran He, Ming-Zhong Ai, Jin-Ming Cui, Yun-Feng Huang, Yong-Jian Han, Chuan-Feng Li, Guang-Can Guo, G. Sierra, C. E. Creffield, “Riemann zeros from Floquet engineering a trapped-ion qubit,” *npj Quantum Information* **7**, 109 (2021), 2021-07-14，DOI [10.1038/s41534-021-00446-7](https://www.nature.com/articles/s41534-021-00446-7)。旧稿的 B.-W. Li / PRL 是占位引文，应删除。

这是把含 zeta 函数的已设计驱动输入离子比特，再以 Floquet 准能交叉和动力学冻结定位零点的实验。它不等价于找到一个自然界固有的 Hilbert–Pólya 哈密顿量。扫描变量 E 是 zeta 的虚部参数，不是零点序号，更不是宇宙年龄。原文称其扫描可进一步延伸，没有把前 80 个零点定义为原理极限。

公开 [Supplementary Information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41534-021-00446-7/MediaObjects/41534_2021_446_MOESM1_ESM.pdf) 的 Tables I、II（PDF 第 3、4 页）可直接转录，已保存并通过 PDF 结构预检。共 269 条估计：序号 1–29 有 Ω=5,8,12,16；序号 30–80 有 Ω=8,12,16。

最后一行的正确解释如下。真值列出版舍入为 201.265。

| Ω/J | 出版原字符串 | 测得虚部 | 括号表示的标准差 | 对舍入真值的残差 |
|---|---|---:|---:|---:|
| 8 | 199.24(576) | 199.24 | 5.76 | −2.025 |
| 12 | 200.14(5) | 200.14 | 0.05 | −1.125 |
| 16 | 200.47(16) | 200.47 | 0.16 | −0.795 |

因此，“误差达到 576”是对末位括号记法的单位误读。它既不是残差 576，也不是标准差 576。其他两个频率在第 80 个位置的括号误差并未同步爆炸。第 80 个零点的真虚部约 201.265，已经在正文所称 E≤200 的扫描端点之外；不能用这一端点外估计辨认全宇宙共同阈值。

表注把误差定义为对 S 数据扰动后进行 4000 次三次插值所得的统计量，并非完备的仪器系统误差预算。主文误差图中的点是估计值减真值，误差条是该插值统计量，两者必须分开。单个 n=74, Ω=8 的 189.50(0) 保留为“舍入为零”标记，绝不可给予无限似然权重。

已输出 `data/processed/he2021_supplement_zeros.csv`，保留原字符串、虚部、正确 sigma、残差、表号和页码。原始逐次离子测量不在公开下载中；主文声明完整 source data 可向通讯作者索取，不能把这份表格转录称作原始 shots 重分析。

## 最近半年确有相关正式发表，但要区分首次公开日期

Shijie Wei, Yue Zhai, Quanfeng Lu, Wentao Yang, Pan Gao, Chao Wei, Junda Song, Franco Nori, Tao Xin, Guilu Long, “The Riemann Hypothesis manifested in dynamical quantum phase transitions,” *Nature Communications* **17**, 8163 (2026)，正式发表 2026-07-01，DOI [10.1038/s41467-026-74935-8](https://www.nature.com/articles/s41467-026-74935-8)。其前身 [arXiv:2511.11199](https://arxiv.org/abs/2511.11199) 已于 2025-11-14 首次公开；不能将它说成过去半年首次出现的实验成果。

真实实验使用五比特核磁共振装置，对前五个零点附近测相干性。论文报告的定位值为 14.12、20.96、25.09、30.44、32.93，正文未给出这五个根估计的独立标准差。文中达到第 10^12 个零点附近的是数值模拟，非实际硬件观测。

这个工作采用特意设计的对数能谱和控制演化。其实验演化变量的量纲转换为 t=epsilon_0 t_ph/hbar；不能把 t 直接解释为宇宙年龄。它增加了可检验的实验平台和资源分析背景，没有测定宇宙时间决定的硬上限。

作者的 [公开仓库](https://github.com/lqf2025/Riemann-data) 和论文列出的 Zenodo DOI [10.5281/zenodo.20590665](https://doi.org/10.5281/zenodo.20590665) 提供数据与代码。Zenodo 直开本次失败，实际下载来自 GitHub 固定 SHA `a8b0b38202d6f330297212f8755c8df04242f9d6`；`git ls-remote` 也核实该 SHA 是 main 的当前引用。保留了仓库 MIT LICENSE 和 README。

下载采用小范围子集，未取约 300 MB 的高零点理论 GLA 文件。作者 README 区分 `data/exp/` 的已处理实验资料、`data/theory/` 的理论资料以及 `data/data*.npz` 的数值生成资料。本项目遵守这个分类：

- `wei2026_processed_tscan.csv`：β=0.3 与 0.5 的 110 个测量点（每组 5 个零点邻域、每个 11 点），保留 x、y、模的均值与标准差。
- `wei2026_processed_beta_scan.csv`：固定第一个零点附近编码时间，18 个 β 扫描点。
- `wei2026_published_zero_estimates.csv`：正文所报五个根估计，明确根误差未发表，不补造误差条。

这些 NPZ 是作者处理后的均值/标准差，当前下载子集没有原始重复实验信号或完整跨点协方差。`t_grid` 是 5×11，`t_flat` 对应 `t_grid.T.reshape(-1)`；转录脚本已经验证此排序，避免把零点邻域和实验点错配。保存的上游 `Fig4.py` 只供审查，未执行。

## 对修订稿的直接建议

保留“有限实验资源只能提供有限精度、可达高度依协议而变”的研究问题。删除“80 被证明是物理硬上限”和“各频率误差在 80 同时指数发散”的结论。把宇宙年龄对应能力增长作为新增假设，而不是上述两项实验的发现。只在明确控制协议、标定、噪声模型、阈值、样本量和能量/时间预算后，比较可达高度。

建议真实实验图并列显示 2021 的中心残差与插值误差，勿混用；2026 可展示处理后相干性数据，单独表明它没有提供宇宙年龄变量。2021 和 2026 的平台、可观测量及协议不同，不应将“80→5”或“5→10^12”拟合成时代演化曲线。

## 可复现文件与检索范围

`python code/extract_public_experiments.py` 从已下载 PDF 与 NPZ 重建 CSV，并输出 `reports/experimental_data_extraction.json`。PDF 页数、结构预检结果见对应 preflight JSON；下载 URL、字节数、SHA256 见 `data/raw/experimental_sources_download_manifest.json`。

本次检索包括 Riemann zeros / trapped ion / quantum simulator / experiment / 2026、精确论文标题、2026-03-20 至 2026-09-20 日期限定，以及 Nature、APS、arXiv 主站核查。找到上述相关正式发表与公开资料；未验证到比它更新的、可用于本课题的实际高零点硬件测量数据。搜索结果还包含个人仓库、数学数值实验及非同行评审声称，未拿来充当新观测。本报告是有记录的针对性检索，不声称穷尽全部文献。
