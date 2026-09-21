#!/usr/bin/env python3
"""Summarize the frozen centered-response sensitivity experiment, irrespective of outcome."""
from pathlib import Path
import json,numpy as np
from scipy.stats import norm
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/order_response'
s=json.loads((OUT/'model_summary.json').read_text());p=json.loads((OUT/'profile_summary.json').read_text());assert all(len(p[l][k])==13 for l in ['2600','2382'] for k in ['G1','A1'])
f=s['2600']['full'];test=s['2600']['heldout2019_2020'];control=s['2382']['full'];rows=[]
for line,block in s.items():
 for kind,r in block['full'].items():
  rows.append(f"| {line} | {kind} | {r['chi2']:.3f} | {str(round(r['order_offset_m_s'],2))+' ± '+str(round(r['order_offset_sigma_m_s'],2)) if 'order_offset_m_s' in r else '固定为 0'} | {', '.join(r['active_bounds']) or '无'} |")
text=f'''# Fe II 2600 相邻级次差异：残差与零质心非高斯响应检验

本轮针对先前发现的约 −60 m/s 级次差异进行了进一步实际计算。**该差异在两个迹线、各曝光删除控制和本次预先固定的正值、零质心非高斯仪器响应族中仍然存在；现有数据没有把来源明确归因于波长标定或线形。** 改变响应形状没有提供可迁移的明显拟合优势，不能据此宣布测出了真实仪器响应，也不能排除其他常规仪器和模板因素。

## 控制的确定与数据

`results/order_response/model_protocol.json` 在查看本轮模型比较结果前冻结。所有模型使用此前已验证的 17 次 ESPRESSO 原生 S2D 提取光谱、同样的质量/原作者 UPL/大气掩膜以及固定 45 分量气体模板。没有新增残差裁剪，也没有改写任何旧曝光结果。主线为 Fe II 2600，对照为 Fe II 2382。

每次曝光有一个共同速度；每个曝光/级次/迹线独立拟合连续谱振幅、斜率和加性零点。每个级次的仪器宽度及形状在 17 次曝光和两迹线间共享。可选的单个参数 δ 定义为**同一跃迁低索引级次减高索引级次的速度**，不涉及 2374 参考线或不同跃迁的原子物理响应。

这里的“固定气体模板”仍然来自整个合并谱。没有重新拟合全部气体分量，也未引入真实测得的提取响应矩阵。两个级次的掩膜、采样位置与误差权重不同，即使实际波长标定完全相同，共享模板的误差也可能通过不同权重产生不同表观位移；因此这个对比不是只针对波长标定的纯粹实验。所有拟合误差以对角 ERRDATA 和固定模板为条件。改变到本轮的联合共享响应参数后，估计量与此前逐曝光单独拟合再汇总的估计量略有不同，约 −60 m/s 的对照应理解为稳定性，而非独立测量。

## 真正的非高斯形状控制

本轮使用正值 Gaussian 混合核

`K(u) = 0.8 N(u; −0.2 d, σ_g²) + 0.2 N(u; +0.8 d, σ_g²)`，

其中 `d = 1.2 cbrt(a) km/s`、`−1 ≤ a ≤ 1`，`σ_g² = (w/2.354820045)² − 0.16 d²`，宽度参数 `1.6 ≤ w ≤ 2.6 km/s`。

它严格满足归一化、零一阶矩及固定总方差 `(w/2.354820045)²`，三阶中心矩为 `0.165888 a (km/s)³`。因此没有允许核的质心漂移与谱线速度形成完全相同的自由度；形状确实能够变成不对称。`w` 是 Gaussian 等效二阶矩宽度，**非 Gaussian 核的实际半高全宽不一定等于 w**。在这个混合族中，归一化和质心约束不会使谱线位移与不对称形状自动正交，特别是混合、饱和谱线。

四个匹配模型为：G0（Gaussian，无级次位移）、G1（Gaussian，有级次位移）、A0（上述不对称核，无级次位移）、A1（上述不对称核，有级次位移）。A0/A1 各用三组预先指定起点并保存全部尝试。宽度在四种模型中均自由，因而本次形状检验不等同于以前的 Gaussian 宽度检验。

## 实际拟合结果

| 跃迁 | 模型 | χ² | δ（m/s，条件局部标准差） | 最优解碰界参数 |
|---|---|---:|---:|---|
'''+ '\n'.join(rows)+f'''

2600 中，G0→G1 的实际改善为 Δχ²={f['G0']['chi2']-f['G1']['chi2']:.3f}；A0→A1 为 Δχ²={f['A0']['chi2']-f['A1']['chi2']:.3f}。允许不对称后级次偏移仍约 {f['A1']['order_offset_m_s']:.2f} m/s，没有被本族消除。相对于 G1，A1 新增两个形状参数只改善 χ²={f['G1']['chi2']-f['A1']['chi2']:.3f}。

A0 在强制 δ=0 时，一个级次的形状参数碰到边界。这说明无位移拟合已用到所选族的极端形状；不能把“本族解释不了”推成“所有不对称响应都已排除”。A1 中 δ 与两个形状参数的局部相关分别为 {f['A1']['offset_correlations']['asymmetry_52']:.3f} 和 {f['A1']['offset_correlations']['asymmetry_53']:.3f}，仍有形状/速度混淆。2600 的两个形状参数仅为 −0.353±0.516、−0.925±0.775（条件局部近似），并未被精确识别；后者虽然优化器未标记为恰好碰界，其不确定范围仍明显受到边界影响。A1 的绝对名义 χ²/自由度在 2600、2382 分别约为 1.34、1.29，因此四模型之间的比较也不代表任何一个模型已充分描述残差及噪声。

在查看主拟合结果后，另计算了 δ 从 −120 到 +120 m/s 的 13 点条件目标函数曲线；每一点重新优化其余所有自由参数。它用于显示实际非线性曲率和边界，并非从局部 Gaussian 误差直接推出新物理显著性。混合参数在 a=0 附近一般不具有连续二阶导数、部分解又碰界，因此本报告不给 Wilks 发现概率或“σ 发现”标签。

2382 对照的级次偏移在 Gaussian 与非 Gaussian 模型下分别为 {control['G1']['order_offset_m_s']:.2f}±{control['G1']['order_offset_sigma_m_s']:.2f} 和 {control['A1']['order_offset_m_s']:.2f}±{control['A1']['order_offset_sigma_m_s']:.2f} m/s，未表现出类似 2600 的约 60 m/s 差异。

## 条件响应迁移检验

只用 2018 年 9 次曝光拟合仪器宽度、形状和可选级次偏移，然后将这些量固定，在 2019–2020 年 8 次曝光中只重新拟合共同曝光速度及逐行连续谱。2600 的测试集 χ²：

| 模型 | 条件测试集 χ² | 相对 G0 |
|---|---:|---:|
'''+ '\n'.join(f"| {k} | {r['chi2']:.3f} | {r['chi2']-test['G0']['chi2']:+.3f} |" for k,r in test.items())+f'''

A1 在整套数据上能降低训练目标函数，但这个响应迁移检验比 G1 **差 {test['A1']['chi2']-test['G1']['chi2']:.3f} χ²**，没有显示新增形状参数的明确可迁移优势。G1 相对 G0 在测试集改善 {test['G0']['chi2']-test['G1']['chi2']:.3f} χ²，规模较小。

此处测试集并非完全未见数据：固定气体模板来自全体曝光合并谱，包含后期曝光的信息。因此这是**在固定共享气体模板条件下的响应迁移诊断**，不是完全独立的预测验证或天体物理确认。真实仪器也可能跨观测期变化，单一共享形状近似不能排除此可能。

## 残差和曝光分组提供的线索

并行的预定义残差诊断保存在 `reports/order_response_residuals_cn.md`。它检查了日期、BERV、信噪比、探测器像素位置/相位，以及分别拟合的两条迹线。主要可复核现象是：

- 2600 的低减高级次差在两迹线中均约 −60 m/s；没有被单个曝光或单个迹线主导。
- 删除任意一次曝光，级次差仍为负；2018 对后期的变化及与 BERV 的斜率都未显示明确偏离零。
- BERV 和谱线在探测器上的位置高度共线，且与观测批次强相关；简单回归不能区分波长位置响应、观测季节或校准批次。
- 强吸收核心的残差幅度大于连续谱，表明固定气体形状/误差模型仍有不足；两级次共有的谱线残差不能直接被认作某一个级次的标定畸变。

## 共同有效覆盖控制

在看到原始级次结果后，又固定了一项只依赖波长几何及原有质量掩膜的控制：每个原生像素的整个波长箱，必须被另一光栅级次原先有效波长箱的并集完全覆盖；不依赖本次拟合残差或通量值选择。采用一次交集，不反复迭代剥除网格边缘。68 个固定模板重拟合全部正常终止且不碰界。

该控制保留 62,154/62,404 个两条线样本，仅去掉 2600 的 63 个、2382 的 187 个像素。2600 的级次差由 −60.247±16.720 变为 −58.650±16.745 m/s；2382 由 −10.528±21.869 变为 −15.445±21.980 m/s。由于两套估计共享绝大部分光子，这些前后变化只作描述，不给独立差异显著性。

在这一具体控制下，缺失覆盖差别没有消除 2600 偏移。但它没有强制两个级次采用相同的采样、信噪权重或响应函数，边缘原生波长箱也不完全相同，因此不能排除通过这些因素表现的模板误差。协议、所选掩膜、拟合及汇总位于 `results/order_response/residual_common_support_*.json`。

## 可以与不可以得出的结论

本次得到了三点更具体的判断：第一，先前级次差并非明显由单次曝光、单迹线或单一观测期造成；第二，对称宽度自由度之外，本次正值零质心形状族也没有消除差异；第三，新增形状参数没有在条件后期迁移检验中胜过 Gaussian 加级次偏移模型。

这些结果**没有确定是波长标定偏差，也没有测出真实非 Gaussian LSF**。可能剩余的原因包括未覆盖的真实响应形状、提取/误差协方差、探测器位置效应以及共享气体模板近似。要区分这些，需要与 2018–2020 实际曝光相匹配的校准响应或独立提取重处理。该吸收体仍只有一个宇宙红移，不提供对数时间律的识别。

## 独立审计

模型与响应族审计 7,278/7,278 项通过，覆盖全部 40 个主拟合和 52 个固定级次偏移目标函数点，包括连续核矩、直接卷积、独立最小二乘复算、解析/有限差分导数、源文件绑定、嵌套模型、条件响应迁移及实际数值加密重拟合。加密网格和提高到 21 点求积后，两条线的最优 δ 改变均小于 0.008 m/s。

残差与几何控制另有 14,533/14,533 项独立检查通过，重建了 136 个迹线拟合和 68 个共同支持拟合；37 个独立标量重优化与存档位移最大相差 0.001022 m/s。另一个只读审查独立重建共同支持的所有掩膜，982/982 项通过。

这些检查证明实现和可复算性，没有验证真实仪器响应、未知像素协方差、模板独立性或新物理。

## 复算与文件

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/order_response_model.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/order_response_profiles.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/order_response_residuals.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/order_response_residuals.py --common-support
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/order_response_validate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/order_response_residual_validate.py
python code/order_response_report.py
```

- `model_protocol.json`：在新比较输出前冻结的模型与控制。
- `model_summary.json`：全数据、2018 训练和后期响应迁移检验的选择结果；全部起点及拟合另有独立 JSON/NPZ。
- `profile_summary.json`：52 个条件 δ 曲线拟合。
- `order_response_models.pdf/.png`：条件目标函数、形状族示意和后期迁移对照。
- `reports/order_response_residuals_cn.md`：独立残差与观察条件诊断。
- `reports/order_response_validation.md`、`results/order_response/validation.json`：模型/核/数值独立验证，以最终审计状态为准。
'''
(ROOT/'reports/order_response_results_cn.md').write_text(text)
fig,axs=plt.subplots(2,2,figsize=(12,8))
colors={'G1':'#286c9b','A1':'#ad582e'}
for ax,line in [(axs[0,0],'2600'),(axs[1,0],'2382')]:
 for kind in ['G1','A1']:
  q=p[line][kind];ax.plot([x['order_offset_m_s'] for x in q],[x['delta_chi2'] for x in q],'o-',ms=3,label='Gaussian' if kind=='G1' else 'Centered asymmetric',color=colors[kind])
 ax.axvline(0,c='gray',lw=.8);ax.set(xlabel='Lower − upper order offset [m/s]',ylabel='Profile objective above own minimum',title=f'Fe II {line}: conditional offset profiles');ax.legend(fontsize=9);ax.grid(alpha=.15)
v=np.linspace(-3.5,3.5,1200);ax=axs[0,1];pars=f['A1']['parameters']
for i,order in enumerate([52,53]):
 w=pars[f'width_{order}'];a=pars[f'asymmetry_{order}'];d=1.2*np.cbrt(a);sg=np.sqrt((w/2.354820045)**2-.16*d*d);kernel=.8*norm.pdf(v,-.2*d,sg)+.2*norm.pdf(v,.8*d,sg);ax.plot(v,kernel,label=f'Order index {order}',color=['#365c86','#985132'][i]);ax.plot(v,norm.pdf(v,0,w/2.354820045),ls=':',alpha=.55,color=['#365c86','#985132'][i])
ax.axvline(0,c='gray',lw=.7);ax.set(xlabel='Kernel velocity [km/s]',ylabel='Normalized response density',title='Assumed family fitted to Fe II 2600');ax.legend(fontsize=9);ax.text(.02,.95,'Unit area; centroid = 0\nDotted: same-RMS Gaussian',transform=ax.transAxes,va='top',fontsize=9);ax.grid(alpha=.15)
ax=axs[1,1];kinds=['G0','G1','A0','A1'];delta=[test[k]['chi2']-test['G0']['chi2'] for k in kinds];ax.bar(kinds,delta,color=['#9da7af','#286c9b','#c5a18a','#ad582e']);ax.axhline(0,color='black',lw=.8);ax.set(ylabel='Later-exposure objective − G0',title='2018 response transferred to 2019–2020');ax.grid(axis='y',alpha=.15)
fig.suptitle('Testing the Fe II 2600 order discrepancy',fontsize=15);fig.text(.5,.012,'Gas template uses all exposures: conditional response transfer, not independent validation. Lower objective is better.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.04,1,.96]);fig.savefig(OUT/'order_response_models.png',dpi=180);fig.savefig(OUT/'order_response_models.pdf');plt.close(fig)
