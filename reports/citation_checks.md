# 引文核查记录

核查日期：2026-09-20。对应文献库：`paper/references.bib`。

本记录区分出版元数据核验、实际可读取的正文范围及论文可以支持的论点。出版商页面、作者预印本和 Crossref 出版元数据的可访问性不同，不能把核验题名及 DOI 写成已经逐页复核全文。本轮没有进行完整撤稿数据库检索，因此不提供“全部文献均无撤稿或更正”的认证。

## 实验与作者已发表论文

| BibTeX 键 | 核查来源与元数据 | 本次证据范围及引用限制 |
|---|---|---|
| `He2021` | [Nature 原文](https://www.nature.com/articles/s41534-021-00446-7)：Ran He 等，*Riemann zeros from Floquet engineering a trapped-ion qubit*，npj Quantum Information **7**, 109 (2021)，在线日期 2021-07-14。 | 从项目已下载的出版商 HTML `citation_*` 标签核对全部作者、题名、卷、文章号、DOI 和日期。最新一次网页工具重开失败，未将失败重开记为成功访问；实验正文及补充表的读取、提取由本项目实验来源报告另记。必须将受驱动量子比特中构造函数的零点估计，与数学零点本身和独立 Hilbert–Pólya 算符构造区分。 |
| `Wei2026` | [Nature 原文](https://www.nature.com/articles/s41467-026-74935-8)：Shijie Wei 等，*The Riemann Hypothesis manifested in dynamical quantum phase transitions*，Nature Communications **17**, 8163 (2026)，在线日期 **2026-07-01**。 | 出版商正文页面可读，缓存 HTML 元数据核对了全部作者和文章信息。[作者预印本](https://arxiv.org/abs/2511.11199)首次提交 **2025-11-14**，预印本题名使用 *Emerges in*，期刊版使用 *manifested in*。因此可称“最近半年正式发表”，不可称“最近半年首次报告的实验”。实际核自旋实验与高阶零点的数值模拟应按原文分开叙述。 |
| `Wang2026Spectral` | [DOI](https://doi.org/10.3390/mca31050193)，用户本地出版 PDF `riemann_hubble/docs/mca-31-00193.pdf`。Liang Wang，*Numerical Spectral Correspondence Between a Non-Autonomous Quadratic Map and the Riemann Zeros: An Exploratory Study*，Mathematical and Computational Applications **31**, 193 (2026)。 | 本轮直接读取 PDF 前两页；第一页明确给出 **2026-09-17** 出版日期、作者、题名、卷及文章号和 DOI。卷期等基础条目同时沿用此前项目文献库；MDPI 网页本轮无法读取，不能据此声称重新核验了在线登记状态。摘要明确把对数冷却称为现象学设定，说明样本外误差较大、模型局部间距不服从 GUE，并否认这是证明或 Hilbert–Pólya 算符构造。新版应引用这些边界，不能仅引用标题来宣称严格同构或物理时间定律。 |

## 数学背景

| BibTeX 键 | 核查来源与元数据 | 引用范围 |
|---|---|---|
| `Montgomery1973` | [出版社 DOI](https://doi.org/10.1090/pspum/024/9944)；本轮 Crossref 该 DOI 返回 Montgomery、题名、1973 年及 181–193 页；卷 24、编辑及书名沿用此前已整理条目。 | 本轮未重新读取原章全文。用于配对相关的历史背景，不能写成所有间距统计已被证明与 GUE 完全相同。 |
| `Odlyzko1987` | [出版社 DOI](https://doi.org/10.1090/S0025-5718-1987-0866115-0)；[作者个人论文列表](https://www-users.cse.umn.edu/~odlyzko/doc/old/zeta.html)列出 Math. Comp. **48** (1987), 273–308。 | 题名、年份、期刊卷页与此前已核对文献库一致。出版商重开受限；本轮不声称逐页重读。用于零点间距的数值检验，不是物理系统的宇宙时间漂移证据。 |
| `BerryKeating1999` | [出版社 DOI](https://doi.org/10.1137/S0036144598347497)；本轮 Crossref 核对作者、题名、SIAM Review **41**(2), 236–266 (1999)。 | 作为谱对应与特征值渐近的研究背景；本轮元数据核验，不是对全文论证的新复核。 |
| `PlattTrudgian2021` | [作者预印本](https://arxiv.org/abs/2004.09765)，[期刊 DOI](https://doi.org/10.1112/blms.12460)；Crossref 核对 Bulletin of the London Mathematical Society **53**(3), 792–797 (2021)。 | 作者摘要明确结果为采用区间算术对 **高度** $3\times10^{12}$ 以下零点的严格数值验证。这里的 $3\times10^{12}$ 是虚部上界，不能误写为零点序号，也不是模拟量子硬件直接观测到了该高度。预印本首次提交为 2020-04-21。 |
| `NISTDLMF` | [NIST DLMF §25.10](https://dlmf.nist.gov/25.10)，在线页面显示版本 1.2.8，发布日期 2026-09-15。 | 本轮读取了零点分布、Hardy $Z$ 函数、计数思想与 Riemann–Siegel 公式；该节明确 $m=\lfloor\sqrt{t/(2\pi)}\rfloor$，此 $m$ 是和式截断长度，不能与零点个数混用。该页没有直接展示新版拟用的完整 Riemann–von Mangoldt 渐近式，不能把它标成该完整公式的精确网页位置。 |
| `Titchmarsh1986` | [NIST 官方书目](https://dlmf.nist.gov/bib/T)核对 *The Theory of the Riemann Zeta-Function*，第二版，1986，D. R. Heath-Brown 编辑，ISBN 0-19-853369-1。 | 补充为标准解析数论背景书目。本轮只核对书目信息，没有重新读取整本书或指定定理页；正文如需要精确页码，应在取得对应页后添加。 |

## 测量与计算资源

| BibTeX 键 | 核查来源与元数据 | 引用范围 |
|---|---|---|
| `BraunsteinCaves1994` | [APS 页面](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.72.3439)：Samuel L. Braunstein、Carlton M. Caves，Phys. Rev. Lett. **72**, 3439 (1994)，1994-05-30。 | 出版元数据及摘要可读，全文需要订阅。引用量子态统计可区分性与测量几何；新版自行给出的 Cramér–Rao 推导须保留局部无偏等条件，不能把统计下界当作所有协议都能达到的精度。 |
| `Giovannetti2006` | [APS 页面](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.010401)：Vittorio Giovannetti、Seth Lloyd、Lorenzo Maccone，*Quantum Metrology*，Phys. Rev. Lett. **96**, 010401 (2006)。 | 元数据及摘要可读，全文访问受限。支持把量子估计精度与资源和策略关联；不直接给出黎曼零点的普适高度上限。 |
| `Bekenstein1981` | [APS 页面](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.23.287)：Jacob D. Bekenstein，Phys. Rev. D **23**, 287 (1981)。 | 元数据及摘要可读。涉及有给定尺寸和能量系统的熵界；适用物理前提必须保留，不能未经论证转为最大素数、最大零点序号或单次实验噪声。 |
| `MargolusLevitin1998` | [作者已发表版本预印本页面](https://arxiv.org/abs/quant-ph/9710043)及 Crossref 核对：Norman Margolus、Lev B. Levitin，Physica D **120**(1–2), 188–195 (1998)。 | **纠正输入中的 DOI**：正确值为 `10.1016/S0167-2789(98)00054-2`，不是 `10.1016/S0167-2789(97)00454-2`。摘要表述的是孤立系统动力学演化速率与高于基态的平均能量的关系，不是把宇宙年龄直接除以普朗克时间就得到可测零点个数。 |
| `Lloyd2002` | [APS 页面](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.88.237901)：Seth Lloyd，Phys. Rev. Lett. **88**, 237901 (2002)，2002-05-24；Crossref 复核。 | 元数据和摘要可读。信息量和操作次数估计可作为资源约束背景，但不提供本项目旧稿的 $T\simeq4200$ 阈值或其时间演化规律。 |
| `Gefen2019` | [Nature 原文](https://www.nature.com/articles/s41467-019-12817-y)、[作者预印本](https://arxiv.org/abs/1811.01762)：Tuvia Gefen、Amit Rotem、Alex Retzker，Nature Communications **10**, 4992 (2019)。出版商 HTML 元数据和 Crossref 一致。 | 网页工具一度跳转失败，直接读取出版商 HTML 元数据成功；作者摘要可读。支持区分特定傅里叶/Rayleigh 分辨准则与更一般量子估计协议，不等于零资源无限精度。 |
| `Boss2017` | [作者预印本](https://arxiv.org/abs/1706.01754)、[Science DOI](https://doi.org/10.1126/science.aam7009)，Crossref 核对 J. M. Boss、K. S. Cujia、J. Zopes、C. L. Degen，Science **356**(6340), 837–840 (2017)。 | 作者摘要说明连续采样的频率分辨率可超越单个量子比特探针的相干时间限制，但依赖外部同步时钟稳定性。因此可以反驳“单次相干时间自动构成一切协议的分辨上限”，不可写成绝对无条件的任意精度。 |

## 本次核查不支持的推论

- GUE 型局部统计相似，不等于零点与所有重原子或原子核能级逐一相同，也不自动产生宇宙年龄依赖。
- 已发表的非自治映射中的迭代变量，尚无实验校准能把它直接等同于宇宙时间。
- 有限资源可以限制指定实验的估计能力；这些引文没有证明数学上的黎曼零点序列在某一高度终止。
- 2026 年新发表实验增加了可分析的平台与数据，并未直接测量宇宙变老导致的零点可观测高度变化。

引用完整性以正文实际引用的键及最终 BibTeX 编译为准；此文件记录来源核查，不代替数值复现和全文同行审阅。
