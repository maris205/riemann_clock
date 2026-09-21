# HARPS HE0515−4414 公共数据可获取性核查

核查时间：2026-09-21 UTC。用途：为 ESPRESSO 相对谱线位移的跨仪器复核寻找独立 HARPS 观测。本次未下载到可验证的 HARPS 科学光谱数组，因此没有 HARPS 重拟合、位移测量或独立确认结果。

## 2021 年论文确实声明公开数据

Milaković et al. (2021), *A new era of fine structure constant measurements at high redshift*, MNRAS 500, 1–21，DOI [10.1093/mnras/staa3217](https://doi.org/10.1093/mnras/staa3217)。原始论文的 Data Availability 明确说明在线补充材料含 LFC 和 ThAr 标定的合并光谱、最终拟合模型及 VPFIT / AI-VPFIT 输入和原子数据。因此不能把本次无法取得附件表述成“作者未公开”。

可读的期刊正文镜像是 [Oxford CDN article-minimal](https://oup.silverchair-cdn.com/article-minimal/5935255)。本地全文保存为 `data/raw/cross_instrument/harps_source_search/oup_minimal.txt`。原始正常期刊页面及几个页面变体在当前环境均返回 HTTP 403。期刊 [Volume 500 Issue 1](https://academic.oup.com/mnras/issue/500/1?browseBy=volume) 显示 Supplementary data，但其链接返回同一文章页面，未给出可验证附件 URL。article-minimal 正文中只有补充材料名称 `supplementary_material`，没有附件链接。若干根据常见出版格式形成的附件候选 URL 返回 CDN `MissingKey`，这仅说明没有可用签名，不能据此判定候选文件名存在或者正确。

## 独立的目录检查

- [INAF 官方仓储](https://openaccess.inaf.it/entities/publication/f0d2f88f-3aa5-4d62-8bae-692bdc137aa3)：通过公开 DSpace REST API 检查 ORIGINAL bundle，只有 `staa3217.pdf`，大小 3,075,562 字节；未列合并光谱附件。元数据及四个 bundle 列表已存档。
- [作者网站](https://milakovic.net/) 与 [研究介绍](https://milakovic.net/research.html)：提供 HARPS 仪器线形产品入口，但本次没有找到 HE0515 合并科学光谱下载。
- 作者 [GitHub dmilakovic](https://github.com/dmilakovic) 公共仓库清单及 `harpslfc`、`Submitted`、网站仓库树已经保存。`harpslfc` 主要为处理例程；本次检查的树中未发现 HE0515 合并科学光谱。未运行其代码。
- Zenodo 全文 `HE0515` 检索返回已有的 ESPRESSO 5512490 产品；DataCite 同关键词结果未找到额外 HARPS 合并光谱记录。检索有覆盖范围限制，不能证明别处不存在数据。
- Crossref 元数据给出论文 PDF 链接，但没有附件关系。其响应已保存。

## 2026 年最新处理工作的意义和限制

Milaković et al., [“Perfect” spectra for ESO’s HARPS spectrograph, arXiv:2607.06809v1](https://arxiv.org/html/2607.06809v1)，2026-07-07，确实对同一目标的约 52.5 小时 HARPS 数据展示 HARPERFECT 新提取方法。其分辨率矩阵和独立采样设计与我们当前的仪器线形及像素相关噪声问题直接相关。本次阅读全文 HTML，没有发现可用的科学光谱/分辨率矩阵下载链接或数据发布声明。故可以把它作为近期方法进展引用，暂时不能当作已取得的新观测数据。它重新处理历史观测，也不是 2026 年新增曝光。

## 对当前分析的处理

本次可执行的跨仪器分析仍应使用已有 UVES 与 ESPRESSO 数组，并明确 UVES 旧合并谱的标定和线形限制。不能把 HARPS 论文中总体 Δα/α 与我们两条 Fe II 强线的自由位移直接比较：统计量、谱线集合、分量模型不同。HARPS 后续复核的明确输入应是作者合并谱、误差、遮罩、波长参考系及标定说明，最好再有与数据同版本的真实仪器分辨率信息。

本次只做公开读取和目录检查，没有联系作者、登录账户或发送任何消息；没有声称 HARPS 验证通过。

所有本地搜索材料的 SHA-256 和大小见 `data/raw/cross_instrument/harps_source_search/search_manifest.json`。
