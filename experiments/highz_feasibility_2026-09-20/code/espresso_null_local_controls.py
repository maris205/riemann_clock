#!/usr/bin/env python3
"""Conditional tangent-space score and exact-forward injection controls.

This does not refit a nonlinear model or calibrate a global discovery test. It
projects relative-shift Jacobian columns against the baseline gas/continuum
Jacobian by direct SVD, avoiding a normal-equation pseudoinverse. Monte Carlo
uses independent standard-normal pixel errors, conditional on the fitted null.
"""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import argparse, hashlib, json, time
import numpy as np
from scipy.stats import chi2
from espresso_conventional_null import Model, ROOT, OUT

SEED = 2609201741

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def project_controls(J, r, n_nuisance, relative_cutoff):
    nuisance, shifts = J[:, :n_nuisance], J[:, n_nuisance:]
    U, singular, nuisance_vt = np.linalg.svd(nuisance, full_matrices=False)
    rank_nuisance = int(np.sum(singular > singular[0] * relative_cutoff))
    basis = U[:, :rank_nuisance]
    nuisance_step = -(nuisance_vt[:rank_nuisance].T / singular[:rank_nuisance]) @ (basis.T @ r)
    effective = shifts - basis @ (basis.T @ shifts)
    r_eff = r - basis @ (basis.T @ r)
    Us, ss, VsT = np.linalg.svd(effective, full_matrices=False)
    rank_shift = int(np.sum(ss > ss[0] * relative_cutoff))
    Us = Us[:, :rank_shift]
    inverse = (VsT[:rank_shift].T / ss[:rank_shift]) @ Us.T
    covariance = (VsT[:rank_shift].T / ss[:rank_shift]**2) @ VsT[:rank_shift]
    estimate = -inverse @ r_eff
    score_coordinates = Us.T @ r_eff
    score = float(score_coordinates @ score_coordinates)
    se = np.sqrt(np.diag(covariance))
    result = dict(relative_singular_value_cutoff=relative_cutoff,
        nuisance_rank=rank_nuisance, nuisance_columns=n_nuisance,
        nuisance_singular_values=singular.tolist(),
        effective_shift_rank=rank_shift, effective_shift_columns=shifts.shape[1],
        effective_shift_singular_values=ss.tolist(),
        local_score_delta_chi2=score,
        conditional_chi2_survival=float(chi2.sf(score, rank_shift)),
        local_shift_estimates_km_s=estimate.tolist(),
        local_shift_sigma_km_s=se.tolist(),
        local_conditional_covariance_km2_s2=covariance.tolist(),
        effective_shift_correlation=(covariance / np.outer(se, se)).tolist(),
        max_nuisance_shift_basis_overlap=float(np.max(np.abs(basis.T @ Us))),
        raw_residual_chi2=float(r @ r),
        nuisance_tangent_projected_residual_chi2=float(r_eff @ r_eff),
        interpretation='Linear tangent-space score with diagonal known Gaussian pixel noise; this is not the nonlinear likelihood-ratio result or a discovery significance.')
    return result, dict(basis=basis, shift_basis=Us, inverse=inverse,
                        covariance=covariance, effective=effective, nuisance_step=nuisance_step)

def validate_projection_algebra():
    """Cross-check projection against a joint linear least-squares solve."""
    rng = np.random.default_rng(726308)
    nuisance = rng.normal(size=(127, 8))
    shifts = rng.normal(size=(127, 5)) + nuisance @ rng.normal(size=(8, 5))
    known = np.array([.01, -.02, .03, -.04, .05])
    residual = nuisance @ rng.normal(size=8) - shifts @ known
    joint = np.column_stack([nuisance, shifts])
    result, _ = project_controls(joint, residual, 8, 1e-12)
    full_solution = np.linalg.lstsq(joint, -residual, rcond=1e-12)[0]
    nuisance_solution = np.linalg.lstsq(nuisance, -residual, rcond=1e-12)[0]
    null_residual = residual + nuisance @ nuisance_solution
    full_residual = residual + joint @ full_solution
    expected_delta = float(null_residual @ null_residual - full_residual @ full_residual)
    return dict(
        recovered_known_shift_max_abs_error=float(np.max(np.abs(np.array(result['local_shift_estimates_km_s'])-known))),
        joint_lstsq_shift_max_abs_discrepancy=float(np.max(np.abs(np.array(result['local_shift_estimates_km_s'])-full_solution[8:]))),
        joint_lstsq_delta_chi2_abs_discrepancy=abs(result['local_score_delta_chi2']-expected_delta))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', default='null')
    parser.add_argument('--draws', type=int, default=20000)
    parser.add_argument('--output-name', default='local_controls')
    args = parser.parse_args()
    began = time.time()
    source = OUT / f'{args.baseline}.json'
    baseline_source_content = source.read_bytes()
    baseline = json.loads(baseline_source_content)
    baseline_snapshot = OUT/f'{args.output_name}_baseline_snapshot.json'
    baseline_snapshot.write_bytes(baseline_source_content)
    assert not baseline.get('free_shifts', False), 'The baseline must fix all line shifts.'
    checkpoint = 'optimizer_success' not in baseline
    cfg = baseline.get('configuration', {}).copy()
    model = Model(free_shifts=True, **cfg)
    p, lower, upper = model.initial()
    lookup = dict(zip(baseline['labels'], baseline['parameters']))
    for i, name in enumerate(model.labels):
        if name in lookup:
            p[i] = lookup[name]
    n_nuisance = model.ngas + model.ncont
    active_bounds = baseline.get('active_bounds', [model.labels[i] for i in range(n_nuisance) if min(p[i]-lower[i],upper[i]-p[i])<1e-5])
    optimizer_success = baseline.get('optimizer_success')
    optimizer_optimality = baseline.get('optimality')
    optimality_text = '尚未报告' if optimizer_optimality is None else f'{optimizer_optimality:.6g}'
    assert baseline['labels'] == model.labels[:n_nuisance]
    assert model.npar - n_nuisance == 5
    assert np.all(p[n_nuisance:] == 0.)
    r0, J, m0 = model.evaluate(p)
    assert np.isclose(r0 @ r0, baseline['chi2'], atol=1e-6, rtol=1e-10)
    sensitivities = []
    for cutoff in [1e-8, 1e-10, 1e-12]:
        result, arrays = project_controls(J, r0, n_nuisance, cutoff)
        sensitivities.append(result)
        if cutoff == 1e-10:
            central, central_arrays = result, arrays
    # Bound-aware first-order stationarity and the finite size of the
    # unconstrained tangent step: a tiny singular direction can project residual
    # power while requiring an impossible or nonlocal parameter displacement.
    nuisance_J = J[:, :n_nuisance]
    gradient = nuisance_J.T @ r0  # derivative of chi2 / 2
    column_norm = np.linalg.norm(nuisance_J, axis=0)
    at_lower = p[:n_nuisance]-lower[:n_nuisance] < 1e-5
    at_upper = upper[:n_nuisance]-p[:n_nuisance] < 1e-5
    projected_gradient = gradient.copy()
    projected_gradient[at_lower & (gradient>0)] = 0.
    projected_gradient[at_upper & (gradient<0)] = 0.
    scaled_gradient = np.divide(projected_gradient,column_norm,
                                out=np.zeros_like(gradient),where=column_norm>0)
    step = central_arrays['nuisance_step']
    proposed = p[:n_nuisance] + step
    violate = (proposed<lower[:n_nuisance]) | (proposed>upper[:n_nuisance])
    top_indices = np.argsort(np.abs(scaled_gradient))[-10:][::-1]
    stationarity = dict(
        gradient_definition='J_nuisance.T @ r is the gradient of chi2/2. At lower bounds positive gradients are KKT-compatible; at upper bounds negative gradients are KKT-compatible. Bound tolerance is 1e-5 in each parameter native unit.',
        bound_tolerance_native_units=1e-5,
        max_abs_raw_gradient=float(np.max(np.abs(gradient))),
        max_abs_bound_projected_gradient=float(np.max(np.abs(projected_gradient))),
        max_abs_column_normalized_bound_projected_gradient=float(np.max(np.abs(scaled_gradient))),
        strongest_column_normalized_feasible_gradients=[dict(parameter=model.labels[i],value=float(scaled_gradient[i])) for i in top_indices],
        nuisance_unconstrained_tangent_chi2_improvement=float(r0@r0-central['nuisance_tangent_projected_residual_chi2']),
        unconstrained_tangent_step_violated_bound_count=int(np.sum(violate)),
        unconstrained_tangent_step_violated_bound_parameters=[model.labels[i] for i in np.flatnonzero(violate)],
        unconstrained_tangent_step_max_fraction_of_bound_width=float(np.max(np.abs(step)/(upper[:n_nuisance]-lower[:n_nuisance]))),
        interpretation='Column-normalized projected gradients diagnose first-order stationarity near the current fit. An unconstrained SVD tangent step may demand large, bound-violating changes; its improvement is not an achievable constrained nonlinear chi2 reduction.')
    # Full-pixel noise is generated in bounded batches. A rho sensitivity model,
    # if configured, makes these white innovations rather than raw pixel draws.
    rng = np.random.default_rng(SEED)
    q = central_arrays['shift_basis']
    rank = q.shape[1]
    samples = []
    for start in range(0, args.draws, 200):
        noise = rng.standard_normal((min(200, args.draws-start), model.ndata))
        samples.append(noise @ q)
    samples = np.vstack(samples)
    score_samples = np.sum(samples**2, axis=1)
    mc_estimates = -samples @ (central_arrays['inverse'] @ q).T
    mc_cov = np.cov(mc_estimates, rowvar=False, ddof=1)
    theoretical_sd = np.sqrt(np.diag(central_arrays['covariance']))
    measured_sd = np.std(mc_estimates, axis=0, ddof=1)
    count = int(np.sum(score_samples >= central['local_score_delta_chi2']))
    probability = float((count+1)/(args.draws+1))
    mc = dict(seed=SEED, draws=args.draws, noise_dimension=model.ndata,
        sample_generation='Independent standard Gaussian values at every retained pixel, projected onto the fixed efficient-shift basis.',
        effective_score_degrees_of_freedom=rank,
        sample_score_mean=float(score_samples.mean()), theoretical_score_mean=rank,
        sample_score_variance=float(score_samples.var(ddof=1)), theoretical_score_variance=2*rank,
        score_quantiles={str(k):float(np.quantile(score_samples, k)) for k in [.5,.9,.95,.99]},
        theoretical_score_quantiles={str(k):float(chi2.ppf(k, rank)) for k in [.5,.9,.95,.99]},
        exceedances=count, observed_score_mc_probability_add_one=probability,
        probability_resolution=float(1/(args.draws+1)),
        probability_note='A finite simulation count is a resolution limit, not proof that the probability equals zero. Only the fixed tangent Gaussian null is simulated.',
        estimated_shift_sd_km_s=measured_sd.tolist(),
        theoretical_shift_sd_km_s=theoretical_sd.tolist(),
        fractional_sd_error=(measured_sd/theoretical_sd-1).tolist(),
        estimated_shift_covariance_km2_s2=mc_cov.tolist())
    # Exact nonlinear flux injections, followed only by the fixed linear local
    # estimator. No observational residual is added to these synthetic data.
    target = 2600
    target_index = model.shift_keys.index(target)
    target_column = n_nuisance + target_index
    injections = []
    for amplitude in [-.05,.05,.10,.20]:
        p_inj = p.copy(); p_inj[target_column] += amplitude
        r_inj, _, m_inj = model.evaluate(p_inj)
        exact_signal = r_inj - r0
        tangent_signal = J[:, target_column] * amplitude
        recovered = central_arrays['inverse'] @ exact_signal
        linear_recovered = central_arrays['inverse'] @ tangent_signal
        truth = np.zeros(len(model.shift_keys)); truth[target_index] = amplitude
        unweighted_errors = np.concatenate([(mi-mb)[line['good']] -
              (J[sl,target_column]*amplitude*line['error'][line['good']])
              for line, sl, mi, mb in zip(model.lines, model.slices, m_inj, m0)])
        # For nonzero rho, J is whitened so the above inverse weighting no longer
        # maps to raw normalized flux. Keep that metric null in that case.
        nonlinear_error = exact_signal - tangent_signal
        injections.append(dict(injected_line=target, reference_line=2374,
            injected_shift_km_s=amplitude,
            local_estimator_recovered_shifts_km_s=recovered.tolist(),
            linear_signal_recovered_shifts_km_s=linear_recovered.tolist(),
            true_shifts_km_s=truth.tolist(),
            nonlinear_local_estimator_bias_km_s=(recovered-truth).tolist(),
            target_recovered_km_s=float(recovered[target_index]),
            target_bias_m_s=float(1000*(recovered[target_index]-amplitude)),
            target_bias_in_conditional_sigma=float((recovered[target_index]-amplitude)/theoretical_sd[target_index]),
            max_abs_nonlinear_minus_tangent_whitened_residual=float(np.max(np.abs(nonlinear_error))),
            rms_nonlinear_minus_tangent_whitened_residual=float(np.sqrt(np.mean(nonlinear_error**2))),
            max_abs_nonlinear_minus_tangent_normalized_flux=None if model.rho else float(np.max(np.abs(unweighted_errors))),
            effective_signal_to_noise=float(np.linalg.norm(q.T @ exact_signal)),
            injected_signal_chi2_noncentrality=float(np.sum((q.T @ exact_signal)**2)),
            local_estimator_noise_sigma_km_s=theoretical_sd.tolist(),
            interpretation='Synthetic flux generated by exact nonlinear conventional forward model; recovery uses only the baseline projected local estimator, not a nonlinear refit.'))
    # Analytic derivative sanity check independent of the larger injections.
    eps = 1e-4
    pp=p.copy();pm=p.copy();pp[target_column]+=eps;pm[target_column]-=eps
    fd = (model.fun(pp)-model.fun(pm))/(2*eps)
    relerr = float(np.linalg.norm(fd-J[:,target_column])/np.linalg.norm(J[:,target_column]))
    checks = dict(projection_algebra=validate_projection_algebra(),baseline_chi2_reproduced=True,
        five_shift_directions_identified=bool(central['effective_shift_rank']==5),
        nuisance_projection_max_overlap_below_1e_8=bool(central['max_nuisance_shift_basis_overlap']<1e-8),
        all_mc_sd_fractional_errors_below_0_04=bool(np.max(np.abs(measured_sd/theoretical_sd-1))<.04),
        exact_linear_signal_recovery_max_abs_error=float(max(np.max(np.abs(np.array(i['linear_signal_recovered_shifts_km_s'])-np.array(i['true_shifts_km_s']))) for i in injections)),
        shift_jacobian_finite_difference_relative_norm_error=relerr,
        shift_jacobian_finite_difference_below_1e_6=bool(relerr<1e-6))
    result = dict(baseline=args.baseline,baseline_file=str(source.relative_to(ROOT)),
        baseline_snapshot_file=str(baseline_snapshot.relative_to(ROOT)),
        baseline_status='incomplete_optimizer_checkpoint' if checkpoint else 'completed_fit',
        baseline_optimizer_success=optimizer_success,baseline_optimizer_optimality=optimizer_optimality,
        baseline_active_bounds=active_bounds,
        observed_pixels=model.ndata, nuisance_parameters=n_nuisance,
        shift_keys=model.shift_keys,reference_line=2374,
        primary_relative_svd_cutoff=1e-10,rank_cutoff_sensitivity=sensitivities,
        primary=central,stationarity=stationarity,monte_carlo=mc,exact_forward_injections=injections,checks=checks,
        warnings=[
            'All estimates are conditional on one fitted conventional gas model and its tangent space.',
            'The retained gas model and transition selection were adopted from the published target; this is a scoped reanalysis, not a blind independent survey.',
            'Pixel covariance is not measured. Independent Gaussian noise and a fixed instrumental profile are conditional assumptions.',
            'Bounds and weak nuisance directions can invalidate an unconstrained tangent approximation. Listed active bounds and rank cutoffs must be reviewed.',
            'A large local score alone cannot distinguish wavelength calibration, omitted gas or blends, instrumental profile error, laboratory systematics, or varying constants.',
            'The simulated nuisance tangent projection does not include model-selection uncertainty or refitting finite-noise nonlinear data.',
            'Exact-forward injections test local sensitivity of the selected algorithm; they are not observed physical shifts.',
            'No cosmic-time trend or inverse-log-squared law is inferred from this one absorber.'],
        configuration=cfg,seconds=time.time()-began,
        provenance_sha256={str(baseline_snapshot.relative_to(ROOT)):hashlib.sha256(baseline_source_content).hexdigest(),
                           'code/espresso_conventional_null.py':sha(Path(__file__).with_name('espresso_conventional_null.py')),
                           'code/espresso_null_local_controls.py':sha(__file__)})
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/f'{args.output_name}.json').write_text(json.dumps(result,indent=2)+'\n')
    # Store compact coordinates only; do not create a large redundant MC file.
    np.savez_compressed(OUT/f'{args.output_name}.npz',
        score_samples=score_samples,shift_noise_samples_km_s=mc_estimates,
        shift_keys=np.array(model.shift_keys),shift_basis=q,
        effective_shift_jacobian=central_arrays['effective'],
        local_conditional_covariance=central_arrays['covariance'])
    report = ROOT/'reports'/f'espresso_null_{args.output_name}_cn.md'
    rows = '\n'.join(f'| {key} | {estimate*1000:+.3f} | {sigma*1000:.3f} |' for key,estimate,sigma in zip(model.shift_keys,central['local_shift_estimates_km_s'],theoretical_sd))
    injrows = '\n'.join(f'| {a["injected_shift_km_s"]*1000:+.0f} | {a["target_recovered_km_s"]*1000:+.3f} | {a["target_bias_m_s"]:+.3f} | {a["max_abs_nonlinear_minus_tangent_whitened_residual"]:.4f} |' for a in injections)
    rankrows = '\n'.join(f'| {a["relative_singular_value_cutoff"]:.0e} | {a["nuisance_rank"]} | {a["effective_shift_rank"]} | {a["local_score_delta_chi2"]:.6f} |' for a in sensitivities)
    report.write_text(f'''# ESPRESSO 常规模型的局部检验与注入控制

{'**临时检查：优化尚未结束，以下结果不得作为最终观测结论。**' if checkpoint else ''}

本文件检验固定常规模型附近的统计灵敏度；它不代替完整非线性模型比较，也不把拟合残差解释成新物理。基线为 `{args.baseline}.json`，保留 {model.ndata} 个像素、{n_nuisance} 个气体和连续谱参数。相对位移以 Fe II 2374 为参考。

## 直接 SVD 的局部检验

对误差加权雅可比矩阵的气体和连续谱列做直接 SVD，再把五个相对位移列投影到其正交补。未使用 $J^T J$ 的伪逆进行该检验。标准差未按残差大小缩放。

| 相对奇异值阈值 | 常规参数局部秩 | 位移局部秩 | 局部 Δχ² |
|---|---:|---:|---:|
{rankrows}

主阈值为 10⁻¹⁰。局部 score Δχ² = {central['local_score_delta_chi2']:.6f}；在固定切空间、独立且已知高斯像素误差下，χ²({rank}) 尾概率为 {central['conditional_chi2_survival']:.6g}。这是条件统计量，不能称为整体发现显著性。

| Fe II 线 | 局部相对位移估计 (m/s) | 条件标准差 (m/s) |
|---|---:|---:|
{rows}

基线优化成功标志为 `{optimizer_success}`，最优性指标 {optimality_text}；触及边界的参数为 `{active_bounds}`。未约束的切空间投影不严格处理边界，弱可识别气体参数也会影响其解释。按雅可比列范数归一化并按边界 KKT 方向投影后，最大一阶梯度为 {stationarity["max_abs_column_normalized_bound_projected_gradient"]:.6g}。无约束的常规参数切线更新将越过 {stationarity["unconstrained_tangent_step_violated_bound_count"]} 个参数边界，其最大步长是相应允许区间宽度的 {stationarity["unconstrained_tangent_step_max_fraction_of_bound_width"]:.6g} 倍。因此切线投影所减去的 χ² = {stationarity["nuisance_unconstrained_tangent_chi2_improvement"]:.6g} 不能解释为可实现的受限非线性改进。

## 高斯噪声模拟

以固定种子 {SEED} 生成 {args.draws:,} 组覆盖全部保留像素的独立标准高斯噪声，然后投影到五维有效位移空间。score 均值 {mc['sample_score_mean']:.5f}（理论 {rank}），方差 {mc['sample_score_variance']:.5f}（理论 {2*rank}）。位移标准差与解析预测的最大相对偏差为 {np.max(np.abs(measured_sd/theoretical_sd-1))*100:.3f}%。超过观测 score 的模拟有 {count} 次；加一估计为 {probability:.6g}，模拟概率分辨率为 {1/(args.draws+1):.6g}。

模拟验证固定切空间中的噪声传播，并没有生成再完整重拟合 147 个气体/连续谱参数的独立非线性实验。真实像素协方差未知，因此模拟不能消除观测噪声模型的不确定性。

## 精确前向模型注入与局部估计恢复

把 Fe II 2600 相对 Fe II 2374 的已知位移注入完整非线性 Voigt、仪器卷积和像素积分模型，其他参数固定。以下用**基线局部投影估计器**恢复，不称为完整非线性拟合恢复；合成信号不叠加实际观测残差。

| 注入 (m/s) | 局部恢复 (m/s) | 局部恢复偏差 (m/s) | 最大非线性－切线差 (像素噪声单位) |
|---:|---:|---:|---:|
{injrows}

小步长中心差分验证位移雅可比：相对范数误差 {relerr:.3g}。纯线性注入的最大恢复误差为 {checks['exact_linear_signal_recovery_max_abs_error']:.3g} km/s。注入结果说明这一算法在所选模型附近能否感知相对位移，并非观测发现。

## 解释边界

该结果只适用于一个吸收系统的一组 Fe II 谱线。波长标定、仪器线扩散函数、未建模速度结构或混合吸收、实验室波长/同位素结构均可能造成相对位移。常规模型的边界、结构选择和未测像素协方差未纳入本模拟。既不据此宣称物理常数改变，也不拟合宇宙时间或 $1/\\ln^2 t$。

机器可读结果：`results/espresso_null/{args.output_name}.json`；模拟和有效雅可比：`results/espresso_null/{args.output_name}.npz`。文件内记录随机种子、数值秩、完整协方差、注入偏差和输入/代码 SHA-256。
''')
    print(json.dumps(dict(primary=central,stationarity=stationarity,monte_carlo=mc,exact_forward_injections=injections,checks=checks,report=str(report)),indent=2))

if __name__ == '__main__':
    main()
