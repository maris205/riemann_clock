# ESPRESSO 历史标定与非高斯仪器线形：可执行性核查

本报告核查的是 HE 0515−4414 的 17 次公开原生曝光能否取得**独立仪器标定**，从而继续追查同一 Fe II 2600 跃迁的相邻级次差异。它不是新的宇宙时间拟合，也没有完成 CCD 重新提取或证明约 60 m/s 差异的原因。

**结论：历史原始标定帧及其关联元数据确实可查，已经找到准确的档案编号；但本轮没有找到与这些科学曝光直接匹配、可立即代入的独立非高斯 LSF 核。默认档案主标定还与原论文采用的 LFC 处理不同，不能不加核对地替换。**

## 1. 已有光谱实际采用什么标定

逐个读取全部 17 个科学 FITS 的头部，均为 `espdr/2.2.3`、`SINGLEHR`、`ESO DET BINX=2`、`BINY=1`。头部明确记录 `wave_cal_source=LFC`，A 光纤主波长文件类别为 `WAVE_MATRIX_LFC_FP_A`。原文件还包含 `WAVEDATA_VAC_BARY` 和像素宽度数组；此前分析使用的原生波长并没有丢失，但它们来自原有处理，不是独立标定测量。完整头部摘录保存在 `data/raw/calibration_readiness/science_exposure_header_extracts.json`。

两个代表曝光的 A 光纤波长主文件都名为 `ESPRESSO_LFC_WAVE_MATRIX_A.fits`，但 `DATAMD5` 分别为：

| 科学原始帧 | 头部记录的波长主文件 DATAMD5 |
|---|---|
| `ESPRE.2018-11-04T07:00:00.643` | `9e3cd8f4b1634624f42a9c8e4094719c` |
| `ESPRE.2020-02-28T01:14:22.246` | `e7524cf800158ca0dbb9da7e31a737c2` |

这些是 FITS 中记录的数据签名，不应误作本地完整文件的 SHA256，也不是档案的唯一文件编号。后续若取得候选主文件，应核对其签名、配方、原始输入和漂移修正链。

原论文第 2.2 节明确说明：2018 年使用夜间开始时取得的 LFC；2019 年单次曝光的 LFC 约早 16 小时；2020 年尝试过 2 月 27 日、29 日和 3 月 3 日的 LFC，最终全部七次科学曝光采用 2 月 27 日的 LFC。不同日期还需对应漂移处理。原作者另做过 ThAr/FP 与 LFC 的整体 α 测量对照；我们的局部、固定模板级次检查不是该整体估计量的直接重测。[Murphy et al. 2022，第 2.2 与 4.3.1 节](https://arxiv.org/html/2112.05819)

## 2. 2024—2026 年方法进展能提供什么

Schmidt 与 Bouchy 的 2024 年论文确实表明 ESPRESSO 的仪器线形具有非高斯、不对称以及随级次、光纤、迹线位置变化的结构。他们采用 **2023 年 1 月 24 日、1HR2x1 模式**的白天标定，使用自建流水线和非参数 LSF，改善了两套波长解的一致性。其数据声明指向 ESO 原始档案；本轮核对的论文正文、数据声明及目标数据仓库没有给出能直接匹配我们 2018—2020 年曝光的 LSF 数组或已固定版本的相应自建流水线。[Schmidt & Bouchy 2024](https://doi.org/10.1093/mnras/stae920)，[可读取全文](https://arxiv.org/html/2404.05283v1)

**同一种观测模式不等于同一时期的响应。** 上述实验在 2022 年标定设备升级之后，不能直接把其展示曲线作为历史科学曝光的独立校正。也不能因为后来的概述提及 LFC 曾停机，就断言 2020 年没有标定：目标论文的具体记录及本次档案查询均证实相关帧存在。

2025 年碘吸收池实验使用 2023 年 5 月的资料，展示了通过科学光路检查波长标定的方法，也指出入纤照明几何的重要性；公开数据项目号为 `60.A-9680(A)`。这是可复用的验证思路，不是我们的历史曝光已经接受过该验证。[Schmidt et al. 2025](https://arxiv.org/html/2504.18485v1)

2026 年 3 月发表的 34 GHz EOM 光梳测试使用 **2024 年 12 月**的实验帧，报告了与级次结构相关的标定差异以及光纤注入问题，继续强调响应建模。它进一步说明高稳定性不自动等于绝对准确性，但不能用于认定我们的 60 m/s 差异来自同一个原因。该文注明项目号 `114.28HD.001`，并描述了其自身数据获取渠道的问题；这种说明不应外推成我们所需历史标定不可获取。[Schmidt et al. 2026](https://arxiv.org/html/2603.17908v1)

## 3. 本轮真正查到的 ESO 数据链

已经实际调用 ESO `tap_obs` 和 DataLink/CalSelector，而不仅是列出网页入口：

- 查询科学帧及 2018、2019、2020 年附近的 LFC 原始帧，保存带字段定义、帧号和 `access_url` 的响应。
- 对 2018-11-04 和 2020-02-28 两个代表科学帧分别取得 `Raw2Master`、`Raw2Raw` 的完整列表和关联树。两个 `Raw2Master` 主记录均标注 `certified=true`、`complete=true`；但这些标记只说明档案关联状态，不能替代科学误差验证。
- 默认 `Raw2Master` 返回的是 `WAVE_MATRIX_THAR_FP_A` / `WAVE_MATRIX_FP_THAR_B`，而原公开科学光谱头部采用 LFC 类别。默认 `Raw2Raw` 也未保证重现原作者手工选择的历史 LFC 日期。
- 找到 2020-02-27 多个 LFC 帧，并核查原始头部。准确的 `SINGLEHR`、2×1 候选包括 `ESPRE.2020-02-27T12:29:36.495`、`ESPRE.2020-02-27T14:21:23.803`（`WAVE,LFC,FP`），以及 `ESPRE.2020-02-27T14:24:30.248`（`WAVE,FP,LFC`）。这些是模式匹配候选，尚未通过原科学头部的主文件签名逐一确认为最终输入。

有限搜索一共取得 55 条 LFC 日期窗口记录，直接读取其中 35 个 `SINGLEHR` 原始头部，确认 17 个 2×1 候选。两个科学帧的关联树显示，2018 年代表帧的 `Raw2Raw` 为完整但未认证；2020 年代表帧自动选取 2 月 28 日的 LFC，与原论文指定的 2 月 27 日不同。这些状态均保留在原始 XML 中。

为追查精确的已处理主文件，又检查了五个模式匹配 LFC 原始帧的 DataLink、两套科学主标定的关联树、公开 TAP 表清单和两个现成主波长文件的头部。五个 LFC DataLink 提供原始帧、CalSelector 服务和夜间日志，没有直接关联已提取的 LFC 光谱或所需 LFC 主矩阵；两个现成主文件的 `DATAMD5` 为 `6337aef5219f317c5235306b400f492f`、`6330f2a141c60224a517272f5e9b3bf5`，也不匹配上表。**本轮没有取得精确匹配的独立主矩阵，因此没有下载替代的 THAR/FP 矩阵冒充重现；这不等于证明其他渠道不存在。**

同一位置的 ESO 已处理科学产品亦可检索，例如 `ADP.2021-04-14T13:40:23.131` 对应首个 2018 年科学原始帧，另有二维级次附件。该产品链不能自动替代论文的自选 LFC 处理；本次仅保存其元数据。默认完整关联包为约 2.24 GB 的主标定组或约 9.72—12.03 GB 的原始标定组，未批量下载。

官方 CalSelector 按仪器标定方案自动关联文件，不能通过该接口任意指定自选的标定类别和数量；因此“自动完整”并不意味着“严格重现作者的 LFC 处理”。[ESO CalSelector 说明](https://archive.eso.org/cms/application_support/calselectorInfo.html)

程序入口为 `https://archive.eso.org/tap_obs/sync`，用 `REQUEST=doQuery`、`LANG=ADQL`、`FORMAT=json`。例如：

```sql
SELECT dp_id, date_obs, ins_mode, dp_cat, dp_type, access_url, datalink_url
FROM dbo.raw
WHERE instrument = 'ESPRESSO'
  AND dp_cat = 'CALIB'
  AND dp_type LIKE '%LFC%'
  AND date_obs >= '2020-02-27T00:00:00'
  AND date_obs < '2020-02-28T00:00:00'
```

科学帧的关联元数据可通过如下形式取得，把末项换为 `calSelector_raw2raw` 可查询原始标定关联：

```text
https://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?ESPRE.2020-02-28T01:14:22.246&eso_download=calSelector_raw2master
```

需要逐项读 `semantics`、`eso_category`、`description`、`access_url`，并检查 `ASSOCIATION_TREE`，不能直接下载整个列表后假定重现成功。[ESO 程序化访问](https://archive.eso.org/cms/eso-data/programmatic-access.html)

## 4. 可立即执行的科学边界

当前能可靠保留的是原生曝光的级次、迹线、批次和模型敏感性检验。下一次真正的独立仪器检验，应先固定匹配历史 LFC 帧、提取配方和漂移链，再从标定帧测量不同级次/迹线上的线形，利用这些标定约束去预测科学谱中的差异；不能先用科学谱的差异调出一个 LSF，再把拟合改善称为独立确认。

具体应优先检查 Fe II 2600 对应的 **5590.810—5595.120 Å** 科学波段，原生 S2D 中从 0 开始计数的 **104/105 与 106/107 行**，即级次索引 52 与 53 的两对迹线。获取标定后还需按相应坐标系和漂移转换对准该范围，不能把真空质心波长直接视作未修正标定帧坐标。所需的新增信息是这些历史位置的独立波长/响应约束。

`ORDER_PROFILE_A/B` 是提取过程中使用的级次空间轮廓产品，不能仅凭名称把它当成沿色散方向的完整 LSF。只调整对称 Gaussian 宽度也不等于已经检验非对称响应。

当前官方 ESPRESSO 流水线版本为 **3.6.0，2026-07-30 发布**；本轮保存并查阅其发布说明和 2026-07-27 版手册，未找到其中明确说明已经集成 2024 年非高斯 LSF 方法的条目。这个结果仅限文档核对，不是对完整源码功能缺失的证明，也不能据此声称升级流水线便会解决该级次差异。[ESO 流水线发布页](https://www.eso.org/sci/software/pipelines/espresso/espresso-pipe-recipes.html)，[3.6.0 手册](https://ftp.eso.org/pub/dfs/pipelines/instruments/espresso/espdr-pipeline-manual-3.6.0.pdf)

本轮为有限范围的来源和可执行性核查，没有安装或运行上游流水线，没有下载多 GB 的原始 CCD 关联包，没有联系作者。所有新材料位于 `data/raw/calibration_readiness/`；原科学 FITS、既有拟合和既有结论保持原样。个别出版社请求返回 HTTP 403 已记录，同时成功读取了 arXiv 的作者全文；不能据这一次失败断言数据不公开。

各来源的准确 URL、HTTP 状态、实际抓取时间和 SHA256 位于 `primary_documents_manifest.json`、`eso_requests_manifest.json`、`eso_supplement_requests_manifest.json`、`eso_historical_product_search_manifest.json`；总的工件校验表为 `source_hash_manifest.json`。实际系统 UTC 抓取时间按原样记录，与论文标示的工作日期分开，不人为改写。
