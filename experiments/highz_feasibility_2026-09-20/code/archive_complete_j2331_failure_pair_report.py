#!/usr/bin/env python3
"""Append explicit bounded diagnostic stops/gradients to the new report."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
s=json.loads((ROOT/'results/archive_complete/J233156-090802/failure_pair/summary.json').read_text())
p=ROOT/'reports/archive_complete_j2331_failure_pair_cn.md';text=p.read_text().split('\n## 保存的一阶诊断与非精密位移')[0]
rows=['','## 保存的一阶诊断与非精密位移','', '| 拟合 | nfev | 停止状态 | Coleman–Li | 归一化梯度映射∞ | 列归一化梯度∞ |', '|---|---:|---|---:|---:|---:|']
for label,r in [('原主分析H0',s['null'])]+[(r['name'],r) for r in s['attempts']]:
 q=r['stationarity'];rows.append(f"| {label} | {r['nfev']} | {r['optimizer_status']} / {r['optimizer_success']} | {q['independent_coleman_li']:.6f} | {q['normalized_gradient_mapping_inf']:.6f} | {q['column_normalized_gradient_inf']:.6f} |")
a=s['alternative'];rows+=['',f"仅供诊断的参数值：1608/1611混合谱区相对2374为{a['shifts_m_s']['1608']:+.6f} m/s，2382谱区相对2374为{a['shifts_m_s']['2382']:+.6f} m/s。二者均非精密常数测量，没有计算或声称可用于年龄曲线的置信区间。",'',f"H1每谱区χ²/像素分别为1608：{a['per_line'][0]['chi2_per_pixel']:.6f}、2374：{a['per_line'][1]['chi2_per_pixel']:.6f}、2382：{a['per_line'][2]['chi2_per_pixel']:.6f}。自由位移在这次有限局部搜索中没有把明显失配修复到原定数值标准；这一事实不能替代完整的总体变化检验。",'', 'H0交叉初值没有把原主分析χ²降低超过0.001，因此未触发额外一次H1交叉重拟合。两个新增结果和原主分析均保存，无覆盖。']
p.write_text(text+'\n'.join(rows)+'\n')
