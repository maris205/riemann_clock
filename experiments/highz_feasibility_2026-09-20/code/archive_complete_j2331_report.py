#!/usr/bin/env python3
"""Render the bounded J233156 campaign without promoting inadequate fits."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/archive_complete/J233156-090802';REP=ROOT/'reports'
def read(name):return json.loads((OUT/name).read_text())
def main():
 s=read('ordinary_model_summary.json');cfg=read('selection.json');p=read(s['selected_null']+'.json');a=np.load(OUT/(p['name']+'.npz'))
 fig,axes=plt.subplots(3,2,figsize=(15,9),gridspec_kw={'width_ratios':[3,1]})
 for row,key in enumerate(cfg['lines']):
  v=a[f'{key}_v'];f=a[f'{key}_flux'];e=a[f'{key}_error'];good=a[f'{key}_good'];m=a[f'{key}_model'];ax=axes[row,0];rx=axes[row,1]
  ax.plot(v,f,color='.45',lw=.65,label='Native coadd');ax.fill_between(v,f-e,f+e,color='.8',alpha=.35);ax.plot(v,m,color='#ad371f',lw=1.2,label='Selected ordinary model');ax.set_ylabel(f'Fe II {key}\nNormalized flux');ax.set_ylim(-.13,1.23);ax.grid(alpha=.15)
  rx.plot(v[good],(f[good]-m[good])/e[good],color='#164d68',lw=.65);rx.axhline(0,color='.5',lw=.6);rx.axhline(3,color='.7',lw=.5,ls='--');rx.axhline(-3,color='.7',lw=.5,ls='--');rx.set_ylabel('Residual / expected error');rx.grid(alpha=.15)
  if row==0:ax.legend(fontsize=8,loc='lower right')
 axes[-1,0].set_xlabel('Velocity relative to catalogue z = 2.143 (km/s)');axes[-1,1].set_xlabel('Velocity (km/s)')
 title=f"J233156−090802: {p['ncomp']} shared components, {p['ndata']} unique native pixels; χ²/ν = {p['chi2_per_ndf']:.3f}"
 fig.suptitle(title+'\n1608 region includes same-gas 1611 opacity; ordinary-model adequacy gate '+('PASS' if s['passes_conditional_adequacy_gate'] else 'NOT PASSED'))
 fig.tight_layout();fig.savefig(OUT/'ordinary_profiles.pdf');fig.savefig(OUT/'ordinary_profiles.png',dpi=150);plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(11,4.4))
 for seed in cfg['seeds']:
  records=[read(f'n{n}_null_s{seed}_os9.json') for n in cfg['component_counts']]
  axes[0].plot([x['ncomp'] for x in records],[x['chi2'] for x in records],'o-',label=f'Start {seed}')
  axes[1].plot([x['ncomp'] for x in records],[x['aicc'] for x in records],'o-',label=f'Start {seed}')
 for ax in axes:ax.set_yscale('log');ax.set_xlabel('Shared gas components');ax.grid(alpha=.2);ax.legend()
 axes[0].set_ylabel('Ordinary-model χ²');axes[1].set_ylabel('Ordinary-model AICc');fig.suptitle('Fixed pixels and finite three-start search; lower-dimensional underfitting is retained')
 fig.tight_layout();fig.savefig(OUT/'gas_complexity.pdf');fig.savefig(OUT/'gas_complexity.png',dpi=150);plt.close(fig)
 lines=['# J233156−090802：有限气体模型搜索结果','',f"固定窗口 −480 至 +500 km/s，三个谱区共 {p['ndata']} 个互不重复的原生像素。1608 谱区显式包含相同气体的邻近 1611 不透明度；2374 与 2382 是目标相对位移对。所有权重均为源 FITS expected-fluctuation 对角误差。",'',
 '| 共享气体分量 | 最佳保存零位移 χ² | 名义自由度 | AICc | 优化器成功停止 |','|---:|---:|---:|---:|---|']
 for x in s['by_count']:lines.append(f"| {x['ncomp']} | {x['chi2']:.6f} | {x['nominal_ndf']} | {x['aicc']:.6f} | {x['optimizer_success']} |")
 lines += ['',f"按预先固定的常规模型 AICc 规则选择 {p['ncomp']} 分量，再以 21 倍像素细分重新优化，得到 χ²={p['chi2']:.9f}，名义自由度 {p['nominal_ndf']}，χ²/ν={p['chi2_per_ndf']:.6f}。优化器成功停止状态为 {p['optimizer_success']}，活动边界为 `{p['active_bounds']}`。",'', '| 谱区 | 原生像素 | χ² | χ² / 像素 |','|---|---:|---:|---:|']
 for x in p['per_line']:lines.append(f"| {x['line']} | {x['ndata']} | {x['chi2']:.6f} | {x['chi2_per_pixel']:.6f} |")
 lines += ['', '**预定常规模型充分性门槛：'+('通过。' if s['passes_conditional_adequacy_gate'] else '未通过。')+'** 门槛要求成功停止、总 χ²/名义自由度不超过 1.5、每谱区 χ²/像素不超过 1.8。这只是条件诊断，不是物理模型有效性的证明。', '']
 if not s['passes_conditional_adequacy_gate']:
  lines+=['因此未自动运行额外相对位移备择，也没有输出可以进入时间曲线的 D=2382−2374 测量。不能把当前残差解释成频率漂移：复杂气体结构、局部极小值、标定、响应或额外混叠仍可参与。']
 else:
  pairpath=OUT/'selected_pair.json'
  if pairpath.exists():
   pair=read('selected_pair.json');lines += [f"随后保存了条件备择；Δχ²={pair['delta_chi2']:.6f}，增加 {pair['extra_shifts']} 个谱区位移。是否双侧成功停止：{pair['both_optimizer_success']}。该数值仅为固定架构的条件比较，需结合独立复核和稀疏 profile 诊断，不代表 α 或宇宙年龄测量。"]
 lines += ['', '多个初值之间的差异和有限计数网格已完整保留。ftol/xtol 成功停止不是全局最优证明；stationarity 一阶量与活动边界保存在各拟合 JSON。没有因残差或位移大小删除像素、重选谱线或无限增加分量。', '', '图形：`results/archive_complete/J233156-090802/ordinary_profiles.pdf` 与 `gas_complexity.pdf`。源选择、物理重叠理由及独立检查见同目录 selection/protocol/preflight 文件和 `archive_complete_j2331_selection_cn.md`。']
 (REP/'archive_complete_j2331_results_cn.md').write_text('\n'.join(lines)+'\n')
 print(p['name'],p['chi2'],p['chi2_per_ndf'],s['passes_conditional_adequacy_gate'])
if __name__=='__main__':main()
