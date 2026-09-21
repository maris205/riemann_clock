#!/usr/bin/env python3
"""Summarize completed endpoints without conflating stop criteria and physics."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/uves_completion'

def write_report(summary):
    """Generate a descriptive report; never convert a finite endpoint into sigma."""
    cases=summary['cases'];byname={r['name']:r for r in cases}
    lines=['# UVES 后续：优化终点与仪器宽度边界核查','',
        '本次仅继续已有 UVES 三条 Fe II 谱线的条件检验，使用归档推荐的 expected-fluctuation 误差。'
        '数据仍是历史曝光形成的 511 个像素；没有新增曝光或宇宙时刻。旧的 cross_instrument 结果与报告保留不动。'
        '曝光日期、SQUAD 共加谱中的记录标定修正及数据来源见 [前期来源核查](cross_instrument_results_cn.md)。','',
        '问题是：在同一个 45 分量气体模型、连续谱、零点和逐线高斯仪器宽度下，'
        '继续优化是否改变原有结论；同时，放宽两种假设共有的宽度边界，能否吸收此前的相对线位移。'
        'H0 固定相对实验室线位，H1 仅增加 Fe II 2382 和 2600 相对 2374 的两个位移。','',
        '## 固定范围与复现','',
        '原始宽度倍率范围为 [0.65, 1.25]，扩大范围为 [0.5, 1.4]；'
        '每一组的 H0 与 H1 拥有完全相同的气体和仪器干扰参数。扩大范围是一项敏感性检验，'
        '并不是独立测出的线展宽函数（LSF）先验。全部拟合使用每原始像素 21 个子采样；'
        '数值审计在保存的最终参数上另算 49 个子采样。','',
        '主拟合的 ftol=1e−10、xtol=1e−11、gtol=1e−6，上限为每例 1800 次残差计算。'
        '采用两项独立报告的局部停止诊断：SciPy 的 Coleman–Li 缩放梯度不超过 1e−6，'
        '且单位盒约束梯度映射不超过 1e−5。仅有 ftol/xtol 成功或目标函数变化很小，'
        '均不写成一阶驻点，更不保证全局最优。起点记录、源代码散列、实际参数和原始像素数组均已保存。','',
        '## 保存的终点','',
        '| 情形 | χ² | χ²/(N−k) | 参数数 | 次数 | 优化器状态 | 缩放梯度 | 盒映射 | 严格驻点 |',
        '|---|---:|---:|---:|---:|---|---:|---:|---|']
    for r in cases:
        status=('上限' if r['optimizer_status']==0 else r['message'].replace('`',''))
        lines.append(f"| {r['name']} | {r['chi2']:.6f} | {r['chi2']/(r['ndata']-r['npar']):.4f} | {r['npar']} | {r['nfev']} | {status} | {r['optimality']:.4g} | {r['box_mapping']:.4g} | {'是' if r['stationary'] else '否'} |")
    if any(r['name']=='wide_h0_cross_bounded' for r in cases):
        lines += ['', '另进行一次范围内预先允许的交叉起点核查：扩大范围 H0 的独立续算终点反而'
            '高于原范围 H0，而后者参数在扩大范围内完全可行，因此以保存的 `bounded_h0` 为起点'
            '再算扩大范围 H0，额外上限 400 次。它仅核查同一模型的优化路径，不增加物理自由度；'
            '结果作为 `wide_h0_cross_bounded` 单独保存，不能把它与原范围 H0 的差额全部归因于 LSF。']
    lines += ['', 'χ²/(N−k) 仅为名义描述量；像素相关性、活跃边界和气体模型选择使其不能直接作为校准后的拟合优度检验。',
        '', '以下按参数是否落在相应可行域内，从全部保存终点选取最低目标值配对。'
        '最好的 H0 来自扩大范围的交叉续算，但其三个宽度倍率均在原范围内，'
        '所以同一个可行参数向量可以用于两个范围；这并不表示它已达到任一范围的驻点。'
        '数值差仅描述这些有限优化终点；'
        '未通过严格驻点检查的配对，不用于计算显著性、p 值或“发现”的置信度。','',
        '| 宽度范围 | 选中 H0 / H1 | χ²(H0)−χ²(H1) | 2382 位移 (m/s) | 2600 位移 (m/s) | 两端均驻点 |',
        '|---|---|---:|---:|---:|---|']
    for p in summary['pairs']:
        n=byname[p['null']];a=byname[p['alternative']]
        bounds='[0.5, 1.4]' if p['wide'] else '[0.65, 1.25]'
        lines.append(f"| {bounds} | {n['name']} / {a['name']} | {p['delta_chi2']:.6f} | {a['shifts_m_s']['2382']:.3f} | {a['shifts_m_s']['2600']:.3f} | {'是' if p['stationary_both'] else '否'} |")
    lines += ['', '## LSF 与优化限制','',
        '| 情形 | 2374 宽度倍率 | 2382 宽度倍率 | 2600 宽度倍率 | 活跃边界数 | 最大盒映射参数 |',
        '|---|---:|---:|---:|---:|---|']
    for r in cases:
        s=r['lsf_scale']
        lines.append(f"| {r['name']} | {s['2374']:.6f} | {s['2382']:.6f} | {s['2600']:.6f} | {len(r['active_bounds'])} | {r['worst_mapping_parameter']} |")
    if len(summary['pairs'])==2:
        first,wide=summary['pairs'];a0=byname[first['alternative']];aw=byname[wide['alternative']]
        n0=byname[first['null']];nw=byname[wide['null']]
        shifts=[aw['shifts_m_s'][k]-a0['shifts_m_s'][k] for k in ['2382','2600']]
        lines += ['',f"扩大 LSF 范围后，选中的 H0 目标值下降 {n0['chi2']-nw['chi2']:.6f}，"
            f"H1 下降 {a0['chi2']-aw['chi2']:.6f}；两个相对位移分别改变 {shifts[0]:.3f} 和 {shifts[1]:.3f} m/s。"]
    if cases and not all(r['stationary'] for r in cases):
        nfailed=sum(not r['stationary'] for r in cases)
        lines += ['', f'**{nfailed}/{len(cases)} 个保存终点未达到预先声明的一阶驻点要求。** 因而本轮可以说明已探索的'
            '对称高斯宽度敏感性，却不能给出稳定全局最小值的证明；目标函数差继续保留优化误差。'
            '低于某个终点 χ² 的常规气体/仪器解释仍可能存在。']
    if cases:
        conds=[r['column_normalized_jacobian_condition'] for r in cases]
        lines += ['', f"局部 Jacobian 经各列二范数归一化后的条件数范围为 {min(conds):.3g}–{max(conds):.3g}。"
            '这说明局部存在很弱的参数组合约束；数值满秩并不代表每个气体分量都被精确识别。'
            '该诊断仅描述所用参数化下的局部导数，不能保证另一个气体模型或全局极值的位置。']
    lines += ['', '## 可以和不可以得出的物理结论','',
        '这些是同一吸收体的独立 UVES 光子，气体分量架构和实验室原子表仍与 ESPRESSO 分析共享；'
        '所有气体参数在 UVES 中重新拟合，但本检验并不是模型也独立的盲复现。'
        '扩大对称高斯宽度范围，只检验一种常规解释。它没有覆盖共加谱的完整像素协方差、'
        '实际非高斯/非对称 LSF、残余波长标定误差，或更灵活的气体结构。','',
        '因此，若相对位移方向仍与此前分析相容，只能称为条件性方向相容；'
        '不能声称“已经排除常规解释”，也不能将本表中的位移直接换算成精细结构常数变化。'
        '同一吸收体不会新增宇宙时间采样，当前阶段仍不检验 1/ln²(t) 规律。','',
        '## 独立数值审计','']
    vp=OUT/'validation.json'
    if vp.exists():
        v=json.loads(vp.read_text());lines.append(f"独立验证状态：**{v['status']}，{v['passed']}/{v['count']} 项通过**。"
            '验证通过是实现与保存结果的核对，不等价于优化已驻点或存在新物理。')
        qs=[c for c in v.get('checks',[]) if '21-to-49 quadrature' in c['name']]
        if qs:
            lines.append(f"在实际保存的 {len(qs)} 个终点上，21→49 子采样的最大 |Δχ²| 为 "
                f"{max(abs(c['chi2_difference']) for c in qs):.6g}，最大单像素变化为 "
                f"{max(c['max_sigma_difference'] for c in qs):.6g}σ。")
        lines.append('旧文件保存与独立 CGS-Voigt/FFT 物理计算的细项见 [验证报告](uves_completion_validation.md)。')
    lines += ['', '复现主拟合：`python code/uves_completion_fit.py --name NEW_UNIQUE_NAME --start PATH_TO_SAVED_ENDPOINT`；'
        'H1 增加 `--free`，扩大 LSF 范围增加 `--wide`。默认误差、采样、容差和上限与本报告一致。'
        '每次使用新名称以保留既有快照。验证：`python code/uves_completion_validate.py`；'
        '汇总与本报告：`python code/uves_completion_report.py`。','',
        '对应的完整参数、数组、散列、优化轨迹和 JSON 汇总位于 `results/uves_completion/`；'
        '优化进度图只展示数值目标的变化，不表示物理显著性。','']
    (ROOT/'reports/uves_completion_results_cn.md').write_text('\n'.join(lines))
def summarize():
    cases=[]
    for f in sorted(OUT.glob('*.json')):
        a=json.loads(f.read_text())
        if not all(k in a for k in ['stationarity_criteria','parameters','source_hashes_at_import']):continue
        m=a['stationarity'];idx=m['most_nonstationary_parameter_index']
        j=np.load(f.with_suffix('.npz'))['jacobian'];norm=np.linalg.norm(j,axis=0)
        normalized_sv=np.linalg.svd(j/np.maximum(norm,1e-300),compute_uv=False)
        cases.append(dict(name=a['name'],wide=a['wide_lsf'],free=a['free_shifts'],chi2=a['chi2'],nfev=a['nfev'],ndata=a['ndata'],
            npar=a['npar'],rank=a['numerical_rank'],optimizer_status=a['optimizer_status'],optimizer_success=a['optimizer_success'],
            column_normalized_jacobian_condition=float(normalized_sv[0]/normalized_sv[-1]),
            message=a['message'],optimality=a['optimality'],box_mapping=m['unit_box_gradient_mapping_inf'],
            worst_mapping_parameter=a['labels'][idx],stationary=a['stationary_by_both_criteria'],
            shifts_m_s={k:v*1000 for k,v in a['shifts_km_s'].items()},lsf_scale=a['lsf_scale'],active_bounds=a['active_bounds'],
            json_sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    pairs=[]
    for wide in [False,True]:
        lower,upper=(.5,1.4) if wide else (.65,1.25)
        # All other bounds are identical. A vector found using wider bounds
        # remains a valid feasible candidate for the narrower domain if every
        # fitted width is inside that domain. Do not transfer stationarity flags.
        feasible=[a for a in cases if all(lower-1e-12<=v<=upper+1e-12 for v in a['lsf_scale'].values())]
        null=[a for a in feasible if not a['free']]
        alt=[a for a in feasible if a['free']]
        if null and alt:
            n=min(null,key=lambda x:x['chi2']);a=min(alt,key=lambda x:x['chi2'])
            pairs.append(dict(wide=wide,null=n['name'],alternative=a['name'],delta_chi2=n['chi2']-a['chi2'],
                stationary_both=n['stationary'] and a['stationary'] and n['wide']==wide and a['wide']==wide,
                selection='Lowest saved objective among parameter vectors feasible in this domain; no stationarity transfer between bounds',
                null_from_other_bounds=n['wide']!=wide,alternative_from_other_bounds=a['wide']!=wide,
                candidates_null=[x['name'] for x in null],candidates_alternative=[x['name'] for x in alt]))
    out=dict(cases=cases,pairs=pairs,scope='Expected-fluctuation error array only; three UVES FeII lines; 45-component architecture; no cosmic-time law or calibrated significance.')
    (OUT/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
    if cases:
        fig,axs=plt.subplots(1,2,figsize=(11,4.3),sharey=True)
        for ax,wide in zip(axs,[False,True]):
            for row in [r for r in cases if r['wide']==wide]:
                a=json.loads((OUT/(row['name']+'.json')).read_text());tr=a['trace']
                x=np.r_[0,[t['evaluation'] for t in tr],a['nfev']]
                vals=np.r_[a['initial_chi2'],[t['chi2'] for t in tr],a['chi2']]
                # Difference from the run's captured start. Subtracting the final
                # value would force an artificial endpoint drop on a log scale.
                y=vals-a['initial_chi2']
                ax.plot(x,y,label=row['name'],lw=1.3)
            ax.set_title('Wider LSF bounds [0.5, 1.4]' if wide else 'Original LSF bounds [0.65, 1.25]')
            ax.set_xlabel('Residual evaluations');ax.grid(alpha=.2);ax.legend(fontsize=7)
        axs[0].set_ylabel(r'$\chi^2$ minus value at the captured start')
        fig.suptitle('UVES objective changes: stopping does not imply stationarity',fontsize=12)
        fig.tight_layout();fig.savefig(OUT/'optimization_progress.png',dpi=180);fig.savefig(OUT/'optimization_progress.pdf');plt.close(fig)
    return out
if __name__=='__main__':
    s=summarize();write_report(s);print(json.dumps(s,indent=2))
