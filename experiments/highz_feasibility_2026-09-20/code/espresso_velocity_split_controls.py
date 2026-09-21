#!/usr/bin/env python3
"""Post hoc disjoint-velocity local shift diagnostic; no nonlinear refitting."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import argparse, json, hashlib
import numpy as np
from scipy.stats import chi2
from espresso_conventional_null import Model, ROOT, OUT, ZREF
from espresso_null_local_controls import project_controls


def digest(data):
    return hashlib.sha256(data).hexdigest()


def difference_control(blue,red,cutoff):
    rb, rr = blue['effective_shift_rank'], red['effective_shift_rank']
    if rb!=5 or rr!=5:
        return dict(identifiable=False,reason='At least one region has fewer than five identifiable relative-shift directions.')
    delta=np.array(blue['local_shift_estimates_km_s'])-np.array(red['local_shift_estimates_km_s'])
    covariance=np.array(blue['local_conditional_covariance_km2_s2'])+np.array(red['local_conditional_covariance_km2_s2'])
    singular,U=np.linalg.eigh((covariance+covariance.T)/2)
    keep=singular>singular[-1]*cutoff
    rank=int(np.sum(keep))
    white=(U[:,keep].T@delta)/np.sqrt(singular[keep])
    statistic=float(white@white)
    sigma=np.sqrt(np.diag(covariance))
    return dict(identifiable=bool(rank==5),covariance_rank=rank,
        difference_definition='blue minus red',difference_km_s=delta.tolist(),
        conditional_difference_sigma_km_s=sigma.tolist(),
        conditional_standardized_difference=(delta/sigma).tolist(),
        conditional_covariance_sum_km2_s2=covariance.tolist(),
        covariance_eigenvalues=singular.tolist(),
        joint_difference_chi2=statistic,nominal_degrees_of_freedom=rank,
        conditional_chi2_survival=float(chi2.sf(statistic,rank)),
        interpretation='Disjoint-pixel conditional covariance sum, holding the fitted baseline and its tangent operators fixed. Correlated calibration/model errors, covariance between pixels, nonlinear refitting and post hoc selection are not calibrated.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',default='null')
    args=parser.parse_args()
    source=OUT/f'{args.baseline}.json';source_bytes=source.read_bytes();base=json.loads(source_bytes)
    snapshot=OUT/'velocity_split_baseline_snapshot.json';snapshot.write_bytes(source_bytes)
    source_hashes={str(snapshot.relative_to(ROOT)):digest(source_bytes)}
    for filename in ['espresso_conventional_null.py','espresso_null_local_controls.py','espresso_velocity_split_controls.py']:
        source_hashes['code/'+filename]=digest(Path(__file__).with_name(filename).read_bytes())
    model=Model(free_shifts=True,**base['configuration'])
    p,lower,upper=model.initial()
    lookup=dict(zip(base['labels'],base['parameters']))
    for i,label in enumerate(model.labels):
        if label in lookup:p[i]=lookup[label]
    nn=model.shift_offset
    assert model.labels[:nn]==base['labels']
    assert len(model.shift_keys)==5
    r,J,_=model.evaluate(p)
    assert np.isclose(r@r,base['chi2'],rtol=1e-10,atol=1e-6)
    velocity=np.concatenate([L['v'][L['good']] for L in model.lines])
    line_ids=np.concatenate([np.full(int(L['good'].sum()),L['key'],dtype=int) for L in model.lines])
    masks={'blue':(velocity>=-70)&(velocity<25),
           'red':(velocity>=25)&(velocity<=140)}
    assert not np.any(masks['blue']&masks['red'])
    results={};primary={}
    for region,mask in masks.items():
        variants=[]
        for cutoff in [1e-8,1e-10,1e-12]:
            result,arrays=project_controls(J[mask],r[mask],nn,cutoff)
            proposed=p[:nn]+arrays['nuisance_step']
            violations=(proposed<lower[:nn])|(proposed>upper[:nn])
            result['unconstrained_nuisance_tangent_step_bound_violation_count']=int(violations.sum())
            result['unconstrained_nuisance_tangent_step_max_fraction_of_bound_width']=float(np.max(np.abs(arrays['nuisance_step'])/(upper[:nn]-lower[:nn])))
            variants.append(result)
            if cutoff==1e-10:primary[region]=result
        results[region]=dict(interval_km_s=[-70,25] if region=='blue' else [25,140],
            interval_convention='left closed, right open' if region=='blue' else 'both ends closed',
            ndata=int(mask.sum()),per_line_pixels={str(k):int(np.sum(mask&(line_ids==k))) for k in model.keys},
            actual_retained_velocity_min=float(velocity[mask].min()),actual_retained_velocity_max=float(velocity[mask].max()),
            rank_cutoff_sensitivity=variants,primary=primary[region])
    differences=[]
    for i,cutoff in enumerate([1e-8,1e-10,1e-12]):
        difference=difference_control(results['blue']['rank_cutoff_sensitivity'][i],results['red']['rank_cutoff_sensitivity'][i],cutoff)
        difference['relative_singular_value_cutoff']=cutoff
        differences.append(difference)
    final=dict(baseline=f'{args.baseline}.json',baseline_chi2=base['chi2'],
        baseline_optimizer_success=base['optimizer_success'],baseline_active_bounds=base['active_bounds'],
        reference_redshift=ZREF,velocity_coordinate='c ln(lambda_obs / ((1+z_ref)*oscillator-strength-weighted isotope reference wavelength))',
        reference_line=2374,shift_keys=model.shift_keys,
        primary_relative_svd_cutoff=1e-10,regions=results,
        difference_cutoff_sensitivity=differences,primary_difference=differences[1],
        checks=dict(full_baseline_chi2_reproduced=True,regions_disjoint=True,
            five_effective_shift_directions_in_all_variants=bool(all(a['effective_shift_rank']==5 for reg in results.values() for a in reg['rank_cutoff_sensitivity'])),
            baseline_parameters_unchanged=True),
        selection_status='Post hoc velocity split requested after inspecting whole-window strong-versus-weak transition offsets.',
        method='Subset rows of the full null-model residual and Jacobian, independently project five relative-shift columns against all nuisance columns in each region by direct SVD. No nonlinear fit, gas-model reselection or parameter-bound enforcement.',
        interpretation='A difference between regions is a profile-specific warning under this conditional model. It cannot establish changing atomic constants or a cosmic-time law.',
        limitations=[
            'The same full-window fitted null fixes both local linearization points; no split-region nonlinear optimization is performed.',
            'Independent pixel errors justify summing the conditional covariance matrices only in the fixed tangent approximation. Actual pixel covariance and shared calibration/model errors are not supplied.',
            'Weak gas directions and active bounds make unconstrained projections only local diagnostics. Bound-violating tangent steps are reported.',
            'The split is post hoc, so nominal probabilities are not globally calibrated significance levels.',
            'No inverse-log-squared time law is tested; both regions belong to the same absorber epoch.'],
        provenance_sha256=source_hashes)
    target=OUT/'velocity_split_local_controls.json';target.write_text(json.dumps(final,indent=2)+'\n')
    table=[]
    d=final['primary_difference']
    for i,key in enumerate(model.shift_keys):
        b,rr=primary['blue'],primary['red']
        table.append(f'| {key} | {1000*b["local_shift_estimates_km_s"][i]:+.2f} ± {1000*b["local_shift_sigma_km_s"][i]:.2f} | {1000*rr["local_shift_estimates_km_s"][i]:+.2f} ± {1000*rr["local_shift_sigma_km_s"][i]:.2f} | {1000*d["difference_km_s"][i]:+.2f} ± {1000*d["conditional_difference_sigma_km_s"][i]:.2f} |')
    ranktable=[]
    for i,cutoff in enumerate([1e-8,1e-10,1e-12]):
        b=results['blue']['rank_cutoff_sensitivity'][i];rr=results['red']['rank_cutoff_sensitivity'][i];dd=differences[i]
        ranktable.append(f'| {cutoff:.0e} | {b["nuisance_rank"]} | {rr["nuisance_rank"]} | {b["local_score_delta_chi2"]:.4f} | {rr["local_score_delta_chi2"]:.4f} | {dd["joint_difference_chi2"]:.4f} |')
    report=ROOT/'reports/espresso_velocity_split_local_controls_cn.md'
    report.write_text(f'''# ESPRESSO 速度区间分割：局部相对位移检查

这是查看整段光谱位移后增加的**事后诊断**。使用 `{args.baseline}.json` 的固定常规气体模型，将误差加权残差和雅可比按实际逐线速度坐标分为蓝侧 −70 ≤ v < 25 km/s 与红侧 25 ≤ v ≤ 140 km/s；二者像素不重叠。速度参考红移为 {ZREF}，相对位移均以 Fe II 2374 为参照。蓝侧保留 {results['blue']['ndata']} 个像素，红侧保留 {results['red']['ndata']} 个像素。

在各区间分别以直接 SVD 去除所有常规参数的局部切空间，然后估计五个相对位移。**没有分别重新进行非线性拟合，也没有在投影时约束参数边界**。标准差保持原始独立像素误差尺度，不按残差缩放。

| Fe II 线 | 蓝侧位移（m/s） | 红侧位移（m/s） | 蓝侧−红侧（m/s） |
|---|---:|---:|---:|
{chr(10).join(table)}

Fe II 2382 和 2600 在两侧均偏负，2586 的偏移较小；红侧强线的负偏移数值更大。但五维区间差并未显示可分辨的不一致，因此这次分段**不能把整段残差归结为仅发生在某一侧的局部谱形问题**，也不证明两侧真实位移相同。

在主奇异值阈值 10⁻¹⁰ 下，蓝侧局部 score Δχ² = {primary['blue']['local_score_delta_chi2']:.5f}，红侧为 {primary['red']['local_score_delta_chi2']:.5f}。把两个不重叠区间的**条件协方差**相加，五维位移差的 χ² = {d['joint_difference_chi2']:.5f}，条件 χ²(5) 尾概率为 {d['conditional_chi2_survival']:.6g}。该概率只描述固定切空间及独立已知高斯像素误差，不包括事后分段、参数边界、模型结构选择和共同标定误差，不能用作整体发现显著性。

| SVD 相对阈值 | 蓝侧常规参数秩 | 红侧常规参数秩 | 蓝侧局部 Δχ² | 红侧局部 Δχ² | 区间差条件 χ² |
|---|---:|---:|---:|---:|---:|
{chr(10).join(ranktable)}

主阈值下，无约束的常规参数切线步在蓝侧越过 {primary['blue']['unconstrained_nuisance_tangent_step_bound_violation_count']} 个参数边界，在红侧越过 {primary['red']['unconstrained_nuisance_tangent_step_bound_violation_count']} 个；最大步长分别为允许区间宽度的 {primary['blue']['unconstrained_nuisance_tangent_step_max_fraction_of_bound_width']:.5g} 和 {primary['red']['unconstrained_nuisance_tangent_step_max_fraction_of_bound_width']:.5g} 倍。这再次表明弱方向的投影不能等同于可实现的受限非线性拟合。

区间之间若表现不同，说明残差可能依赖局部速度结构、饱和度或混合吸收，值得作进一步常规模型控制。即使两个区间方向相同，也不能独自排除共同标定或实验室数据系统误差。此处两个区间属于同一吸收系统、同一宇宙时期，不提供时间演化证据，也不检验 $1/\\ln^2 t$。

机器结果保存在 `results/espresso_null/velocity_split_local_controls.json`；该文件记录逐线像素数、所有截断阈值结果、五维完整协方差、基线快照和代码 SHA-256。原有 `local_controls` 结果保持不变。
''')
    print(json.dumps(dict(primary_blue=primary['blue'],primary_red=primary['red'],difference=d,report=str(report)),indent=2))

if __name__=='__main__':main()
