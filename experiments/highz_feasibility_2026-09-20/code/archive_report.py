#!/usr/bin/env python3
"""Render archive-quality report and standalone, data-driven coverage heatmap."""
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
from archive_screen import ROOT,OUT,LINES

def main():
    screen=json.loads((OUT/'catalogue_screen.json').read_text())
    strict=json.loads((OUT/'strict_catalogue_screen.json').read_text())
    source=json.loads((OUT/'complete_source_manifest.json').read_text())
    quality=json.loads((OUT/'absorber_quality.json').read_text())['results']
    summary=json.loads((OUT/'quality_summary.json').read_text())
    quality=sorted(quality,key=lambda r:(r['z_abs'],r['target']))
    mat=np.zeros((len(quality),len(LINES)));cnr=np.zeros_like(mat)
    for i,r in enumerate(quality):
      for j,l in enumerate(r['line_rows']):
        mat[i,j]=3 if l['detected_strict_candidate'] else 2 if l['detected_primary_candidate'] else 1 if l['metadata_eligible'] else 0
        cnr[i,j]=l['native_CNR'] or 0
    colors=['#e7e9ec','#d89c9b','#eac687','#80b7ae']
    fig,ax=plt.subplots(figsize=(11,15))
    ax.imshow(mat,cmap=ListedColormap(colors),norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4),aspect='auto')
    ax.set_xticks(range(len(LINES)),['Fe II '+k for k in LINES],rotation=45,ha='right')
    ax.set_yticks(range(len(quality)),[f'{r["target"]}  z={r["z_abs"]:.3f}'+('  *' if r['passes_computational_readiness'] else '') for r in quality],fontsize=8)
    for i in range(len(quality)):
      for j in range(len(LINES)):
        if cnr[i,j]>0:ax.text(j,i,str(round(cnr[i,j])),ha='center',va='center',fontsize=7,color='#203039')
    ax.set_xticks(np.arange(-.5,len(LINES),1),minor=True);ax.set_yticks(np.arange(-.5,len(quality),1),minor=True)
    ax.grid(which='minor',color='white',linewidth=1);ax.tick_params(which='minor',bottom=False,left=False)
    ax.set_title('Actual native-pixel audit: 36 catalogue-selected DLA absorbers\nCell numbers: median 1/error per native pixel; * = conservative multiplet-readiness gate',fontsize=12,pad=15)
    labels=['Excluded by initial metadata policy','Metadata candidate; actual detection/quality gate failed','Detected; broad H2O-band warning','Detected outside broad H2O bands']
    fig.legend(handles=[Patch(facecolor=c,label=l) for c,l in zip(colors,labels)],loc='lower center',ncol=2,fontsize=8)
    fig.text(.5,.052,'Colours describe computational readiness, not calibrated frequency measurements.\nEW detection uses diagonal errors; weak tellurics, sky, blends, saturation and LSF remain to be checked.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.08,1,1]);fig.savefig(OUT/'archive_line_readiness.pdf');fig.savefig(OUT/'archive_line_readiness.png',dpi=160);plt.close(fig)
    lines=['# 现有 UVES 档案扩展：实际光谱质量分析', '',
      '## Material Passport','',
      '- 类型：真实公开光谱的计算分析与可复现性验证。',
      '- 数据源：Murphy et al. 的 UVES SQUAD DR1，固定 Git 提交 `a0cdc8e7b99f2b01a45d919988af9d60e6447d19`。',
      '- 分析对象：已知 DLA 目录中的 Fe II 多重谱线；不是全部金属吸收系统普查。',
      '- 阶段：读取原生合并光谱、固定规则质量筛选、等效宽度与表观光学深度；没有拟合宇宙时间函数。',
      f'- 报告生成 UTC：{datetime.now(timezone.utc).isoformat()}。时间戳采用计算环境时钟。','',
      '**结论：现有公开数据确实可用于扩大分析。本轮已经完成目录筛选并实际下载、检查全部通过初筛的 36 个 DLA 吸收系统。最终 6 个系统满足保守的多谱线计算准备条件，但没有任何一个因此自动成为经过标定的精细结构常数或宇宙时间测量。**','',
      '## 数据数量与选择边界','',
      '| 阶段 | 数量 | 含义 |','|---|---:|---|',
      f'| 主目录行数 | {screen["master_rows"]} | 包括没有最终合并谱的行 |',
      f'| 有最终光谱的类星体 | {screen["final_spectra_positive_exposure_count"]} | 不是 467 个精密常数测量 |',
      f'| DLA 视线／吸收系统 | {screen["DLA_sightlines"]}／{screen["DLA_absorbers"]} | 同一视线可能有多个吸收红移 |',
      f'| 通过初步目录条件的吸收系统 | {screen["metadata_eligible_absorbers"]} | 覆盖、代理 CNR、森林及主要大气带条件 |',
      f'| 实际下载的不同视线 | {source["n_unique_sightlines"]} | 全部初筛合格系统；下载失败 0 |',
      f'| 实际 Fe II 谱线窗口 | {summary["actual_windows"]} | 每个系统 9 条跃迁，保留全部排除原因 |',
      f'| 排除论文列出全部宽水汽带后的目录候选 | {strict["eligible_count"]} | 保守条件会舍弃某些实际可用窄区间 |',
      f'| 通过原生像素质量与共同吸收条件 | {summary["computational_readiness_count"]} | 进入气体模型开发的候选，非已验证物理变化 |','',
      f'下载归档文件合计 **{source["total_archive_bytes"]/1e6:.2f} MB**。档案为历史曝光经原作者归一化、合并的高分辨率光谱；不是 2026 年新取得的光子观测。波长、归一化流量、误差、预期波动、像素状态和贡献记录均已读取，未对原生数据重新采样。',
      '',
      '来源：[SQUAD 原作者仓库](https://github.com/MTMurphy77/UVES_SQUAD_DR1)、[数据论文 DOI](https://doi.org/10.1093/mnras/sty2834)、[原论文 arXiv](https://arxiv.org/abs/1810.06136)。光谱门户声明 CC BY-SA（未指定版本），代码／目录仓库声明 CC BY 4.0；两者分别保留。','',
      '初筛曾按第三高的合格谱线 CNR 排出前五个目标，随后根据数据论文补充宽水汽带警告并增加三个候选。用户要求把可行分析做完后，下载范围扩展为初筛全部 36 个吸收系统。这些阶段的规则、时间戳、哈希与选择结果全部保存；扩展发生在任何相对位移测量之前，没有按位移方向或拟合收益挑样本。','',
      '## 固定质量条件及其意义','',
      '目录筛选采用吸收红移，而非类星体发射红移：每条线完整 ±250 km/s 窗口须处于目录波长覆盖内、处于 4000–9000 Å、位于类星体 Lyα 红侧 3000 km/s 之外，并避开指定主要大气带；吸收系统与类星体本身也须相隔至少 3000 km/s。至少三条 f≥0.01 的 Fe II 线须达到目录 CNR 代理值 20／2.5 km/s。目录 CNR 只有五个参考波长，不能替代实际流量误差。','',
      '后续严格标记另排除数据论文列出的 H₂O 宽带 6470–6600、6830–7450、7820–8620、8780–10000 Å，并加 ±30 km/s 的地球公转速度余量。这是波长警告与保守剔除，并未使用逐次曝光的大气透射模板。通过该条件也不保证没有弱大气线、天空残差或无关吸收。','',
      '原生像素要求：有效比例至少 98%；全窗口中位 1/误差至少 15／原生像素；固定 ±100 km/s 核心等效宽度除以对角误差传播值至少 5。至少三条合格强跃迁中，须出现至少三个共同的 5 km/s 网格位置达到逐像素名义 3σ 吸收；插值两侧原生像素均须有效。网格位置并非独立平均箱，不将其合计解释为联合显著性。还要求 2374 或 2586 至少一条相对较弱跃迁有五个中等吸收、未达到低流量阈值的像素。','',
      '低流量阈值为 F≤max(0.1,3σ)。表观光学深度使用 max(F,0.02,3σ) 防止对噪声取对数，计算 Na(v)=3.768×10¹⁴ τa/(fλ)，单位 cm⁻²/(km/s)。进入流量下限的像素仅给饱和诊断，不能恢复真实柱密度。等效宽度按 FITS 真实对数网格边界积分；统计误差暂按对角传播，未证明像素独立或连续谱无误差。','',
      '9 条跃迁的实验室数据统一来自已归档 Murphy/Berengut 2014 原子表的陆地同位素组合波长与 f 值。目录红移常只保留三位小数，速度零点可能偏移几十 km/s；这些窗口不能拿来测量精密线心偏移。强、弱线因饱和而出现不同通量形状或表观质心，也不是频率随时间变化的证据。','',
      '## 六个实际候选','',
      '| 目标 | 吸收红移 | 合格强跃迁 | 相对较弱线 | 窗口边缘额外吸收警告 |','|---|---:|---|---|---|']
    for r in quality:
      if r['passes_computational_readiness']:
        lines.append(f'| {r["target"]} | {r["z_abs"]:.3f} | {", ".join(r["strict_detected_strong_lines"])} | {", ".join(r["resolved_weak_line_candidates"])} | {", ".join(r["edge_absorption_warning_lines"]) or "未触发固定阈值"} |')
    lines += ['',
      '其中 J232128−105122（z=1.629）和 J004131−493611（z=2.248）的主要多线吸收较紧凑，已另做具有共同气体参数的探索性物理轮廓拟合，详见 [物理轮廓结果](archive_profile_results_cn.md)。其他候选包含更复杂的速度结构或窗口边缘吸收。这里的“候选”仅表示存在可建模的信息，不是无污染证明。专门的物理试验与本筛选的质量计数分别报告。','',
      '两个新增系统的主拟合使用光谱提供的统计误差；另保留选定气体结构和全部像素，完成了 expected-fluctuation 误差数组的非线性重拟合，详见 [误差选择敏感性](archive_profile_expected_noise_cn.md)。这项控制不重新选择分量数，也不等于已经估计完整像素协方差。','',
      '严格目录候选中的 J084424+124546（z=1.864）没有通过最终门槛：实际 Fe II 1608 CNR 约 11.4／原生像素，低于预定 15，剩余线又受宽水汽带条件限制。这说明目录代理 CNR 确实不能直接当作谱线测量精度。','',
      '## 标定与物理解释限制','',
      '对全部 32 个目标的 UPL 日志另做检查：没有显式记录非零 VSHT 波长平移／斜率校正。这不等于所有其他标定步骤均不存在；它说明本轮没有建立这些一般档案合并谱具备亚 km/s 相对波长精度的证据。原数据论文讨论的长程／局部波长畸变、谱阶合并、blaze 残差、LSF 与大气污染必须继续处理。不能把这些候选直接接到约 100 m/s 的 ESPRESSO 线索上。','',
      '本轮没有给出 Δα/α，也没有将这些谱线的表观柱密度差、饱和差或形状差解释为随时间变化，更没有拟合 1/ln²(t/t*)。需要先规定假说对应的原子响应，再对多个吸收系统完成可比的物理与仪器分析。','',
      '**范围界限：本轮完成的是“155 个已知 DLA 中，符合预先固定条件的全部候选”的实际数据筛查。467 条类星体视线中尚可能有未列入 DLA 目录的金属吸收系统；不能写成公开数据已全部穷尽。**','',
      '## 独立复核','',
      '独立脚本从原始 CSV、FITS 压缩包和原子数据重新构造结果，未导入生成分析代码：目录 5899 项、原始数据解码与严格筛选 939 项、原生窗口数值与筛选 11529 项，合计 **18367 项机械检查全部通过**。检查数量不代表独立观测数量或物理显著性。复核曾发现使用名义像素间距导致约 2–4 ppm 积分尺度差，已改用真实 FITS 对数网格，并修正了共同吸收网格的描述及跨掩码插值限制；最终六候选计数不变。','',
      '## 可复现文件','',
      '- `results/archive_expansion/catalogue_screen.json`：全部 155 个系统、1395 个目录窗口与初筛排除理由。',
      '- `results/archive_expansion/strict_catalogue_screen.json`：宽水汽带附加筛选。',
      '- `results/archive_expansion/complete_source_manifest.json`：32 份真实光谱来源、日期、大小和 SHA256。',
      '- `results/archive_expansion/actual_quality_policy.json`：最终质量规则（保留首版与审查澄清）。',
      '- `results/archive_expansion/window_quality.csv`：324 个原生像素窗口的数值诊断。',
      '- `results/archive_expansion/absorber_quality.json`、`quality_summary.json`：逐系统结论。',
      '- `results/archive_expansion/all_absorber_profiles.pdf`：全部 36 个系统的九线光谱页。',
      '- `results/archive_expansion/candidate_AOD_diagnostics.pdf`：六个候选的表观柱密度对照。',
      '- `results/archive_expansion/archive_line_readiness.pdf`、`.png`：实际质量热力图。',
      '- `reports/archive_quality_review.md`：独立源文件、算法与物理解读复核。','',
      '重新运行顺序：`archive_screen.py` → `archive_fetch.py` → `archive_strict_screen.py` → `archive_expand_all.py` → `archive_quality.py` → `archive_report.py`。原始结果的冻结规则和检索时间戳用于记录本次历史；重跑会产生新的运行时间戳，不应冒充事前注册的全新分析。','']
    report=ROOT/'reports/archive_expansion_results_cn.md';report.write_text('\n'.join(lines))
    print(report)
if __name__=='__main__':main()
