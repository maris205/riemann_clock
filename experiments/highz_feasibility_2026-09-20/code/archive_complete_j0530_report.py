"""Summarize the bounded, frozen-selection J0530 physical-profile campaign."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/archive_complete/J053007-250329'
def read(name):return json.loads((OUT/name).read_text())
def main():
    ordinary=read('ordinary_model_summary.json');pair=read('selected_pair.json') if (OUT/'selected_pair.json').exists() else None
    h0=read((pair['null'] if pair else ordinary['selected_null'])+'.json');z=np.load(OUT/(h0['name']+'.npz'))
    h1=read(pair['alternative']+'.json') if pair else None;a=np.load(OUT/(h1['name']+'.npz')) if h1 else None
    fig,ax=plt.subplots(4,2,figsize=(13,10),sharex=True,gridspec_kw={'width_ratios':[2,1]})
    for row,k in enumerate(h0['configuration']['lines']):
        v=z[f'{k}_v'];f=z[f'{k}_flux'];e=z[f'{k}_error'];m=z[f'{k}_model'];q=z[f'{k}_good']
        ax[row,0].errorbar(v[q],f[q],e[q],fmt='.',ms=2,color='.5',alpha=.6,label='Native coadd and expected error')
        ax[row,0].plot(v,m,color='#1c548b',lw=1.4,label='Shared gas H0')
        ax[row,1].plot(v[q],(f[q]-m[q])/e[q],'.-',ms=2,lw=.6,color='#1c548b',label='H0')
        if a is not None:
            ma=a[f'{k}_model'];ax[row,0].plot(v,ma,color='#d1791d',ls='--',lw=1,label='Relative shifts H1')
            ax[row,1].plot(v[q],(f[q]-ma[q])/e[q],color='#d1791d',ls='--',lw=.8,label='H1')
        ax[row,0].set_ylabel(f'Fe II {k}\nNormalized flux');ax[row,0].axhline(1,color='.7',lw=.5);ax[row,0].set_ylim(-.07,1.10)
        ax[row,1].axhline(0,color='.7',lw=.5);ax[row,1].axhline(3,color='.8',lw=.5,ls=':');ax[row,1].axhline(-3,color='.8',lw=.5,ls=':');ax[row,1].set_ylabel('(data−model)/error')
        for axis in ax[row]:axis.grid(alpha=.15);axis.set_xlim(-110,130)
    ax[0,0].legend(fontsize=8);ax[0,1].legend(fontsize=8)
    for axis in ax[-1]:axis.set_xlabel('Velocity relative to z = 2.141 (km/s)')
    title=f'J053007−250329: {h0["ncomp"]} shared Voigt components, 480 native pixels\nH0 χ² / nominal dof = {h0["chi2"]:.2f} / {h0["nominal_ndf"]}'
    if pair:title+=f'; conditional Δχ² = {pair["delta_chi2"]:.2f} for {pair["extra_shifts"]} additional shifts'
    else:title+='; conventional-model gate not passed'
    fig.suptitle(title);fig.tight_layout();fig.savefig(OUT/'physical_profiles.png',dpi=170);fig.savefig(OUT/'physical_profiles.pdf');plt.close(fig)
    text=['# J053007−250329：固定选区后的物理拟合结果','',
          '这是一项探索性常规模型检验。模型没有精细结构常数参数，也没有宇宙时间或对数时间规律参数。',
          '', '采用 Fe II 1608、1611、2374、2382 的相同原生像素，共 480 个；共同窗口 −110 至 +130 km/s。误差为原始 FITS expected-fluctuation 对角权重。窗口和跃迁集合在相对位移拟合前固定；未按残差剪切像素。','',
          '|共享分量数|H0 χ²|名义自由度|AICc|成功停止|','|---:|---:|---:|---:|:---|']
    for r in ordinary['by_count']:text.append(f'|{r["ncomp"]}|{r["chi2"]:.6f}|{r["nominal_ndf"]}|{r["aicc"]:.6f}|{r["optimizer_success"]}|')
    text+=['',f'按常规模型 AICc 选择 {h0["ncomp"]} 个共享分量，随后以 21 倍像素积分细分优化。最终所用 H0：χ²={h0["chi2"]:.6f}，名义自由度={h0["nominal_ndf"]}，成功停止={h0["optimizer_success"]}。',
           f'固定的条件充分性门槛通过状态：{ordinary["passes_conditional_adequacy_gate"]}。该门槛只判断此轮残差规模是否允许继续条件比较，不证明模型真实，也不排除混合或波长标定问题。','']
    if pair:
        u=pair['local_pair_uncertainty'];sig=u['rank_sensitivity'][1]['conditional_sigma_m_s'];text += [f'H1 在相同气体分量结构上增加 {pair["extra_shifts"]} 个相对位移：χ²={pair["chi2_alternative"]:.6f}，条件改善 Δχ²={pair["delta_chi2"]:.6f}。两模型均成功停止={pair["both_optimizer_success"]}。',
        f'2382 相对于 2374 的位移为 {u["shift_m_s"]:.3f} m/s；局部线性条件标准差为 {sig:.3f} m/s（SVD cutoff=10⁻¹⁰）。此数值固定了所选气体结构、噪声权重和高斯仪器模型；没有合并标定误差、未知混合及分量选择误差。不得解释为 α 或宇宙时间变化测量。',f'相同条件下，H0 的 AICc={h0["aicc"]:.6f}，H1 的 AICc={h1["aicc"]:.6f}；加入位移后的 AICc 增加 {h1["aicc"]-h0["aicc"]:.6f}。因此小幅 χ² 改善没有抵消这一参数惩罚，当前结果没有为额外位移增加有力证据。','']
        if (OUT/'pair_profile.json').exists():
            pp=read('pair_profile.json');text+=['|固定 2382 位移（m/s）|优化后 χ²|相对自由 H1 的 Δχ²|成功停止|','|---:|---:|---:|:---|']
            for r in pp['grid']:text.append(f'|{r["value_m_s"]:.3f}|{r["chi2"]:.6f}|{r["delta_from_free"]:.6f}|{r["optimizer_success"]}|')
            text+=['','这只是稀疏的重新优化剖面，不是经过验证的置信区间。中央固定点比自由 H1 低约 8×10⁻⁶，仅在数值停止精度范围内；没有出现有实质意义的更低剖面点。','']
    text += ['H0 边界命中：'+', '.join(h0['active_bounds'])+'.']
    if h1:text+=['H1 边界命中：'+', '.join(h1['active_bounds'])+'.']
    text += ['', '16 个分量位于预先固定的搜索网格上限，因此本轮没有证明分量架构已经稳定或唯一。3 个初始值和有限迭代只构成有界搜索，不能证明全局最优。优化器按相对目标值或步长成功停止，并不表示所有一阶导数已接近零。', '固定最终参数把积分细分从 21 提至 55 后，H0 和 H1 的 χ² 变化分别约 2.21×10⁻⁵、2.07×10⁻⁵，Δχ² 变化仅 1.43×10⁻⁶；这验证了当前参数下积分误差很小，不替代重新优化或架构稳定性检验。', '', '主要物理限制：2382 有宽黑核，1608 与 2374 也接近饱和；弱线 1611 有助于约束柱密度，但不能唯一分解窄气体成分。各线的连续谱、零点与高斯仪器宽度均有自由度。源光谱的波长标定残差和像素相关并未由本次对角误差拟合完全描述。','', '独立输入选区核查 38/38 项通过；独立拟合与局部误差复核 353/353 项通过，详见 `archive_complete_j0530_selection_review.md` 与 `archive_complete_j0530_fit_review.md`。这些检查验证数据和计算一致性，不证明模型充分或新的物理效应。', '', '拟合图：`results/archive_complete/J053007-250329/physical_profiles.pdf`。原始选区图、源文件哈希、全部初始值和未收敛记录均保留在同一目标目录。']
    (ROOT/'reports/archive_complete_j0530_results_cn.md').write_text('\n'.join(text)+'\n')
if __name__=='__main__':main()
