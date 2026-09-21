#!/usr/bin/env python3
"""Fixed-point information loss from frozen strong-core masks under one GLS kernel.

There is no new optimization or test of the observed residual here. All three
designs use exactly the earlier null_cross physical parameter vector, including
zero relative shifts. Gas/continuum tangent directions are profiled by SVD.
The illustrative signal was specified before this calculation: both strong
lines shift by +100 m/s; the other three free relative shifts remain zero.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.linalg import qr, solve
from espresso_conventional_null import Model, ROOT
from espresso_empirical_gls import GLSModel
from joint_controls_fit import JointModel

OUT = ROOT / 'results/mask_information'
BASELINE = ROOT / 'results/espresso_null/null_cross.json'
CUTOFFS = (1e-8, 1e-10, 1e-12)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    baseline = json.loads(BASELINE.read_text())
    raw = Model(free_shifts=True)
    point, lower, upper = raw.initial()
    lookup = dict(zip(baseline['labels'], baseline['parameters']))
    point = np.array([lookup.get(key, val) for key, val in zip(raw.labels, point)])
    assert raw.shift_offset == len(baseline['parameters']) == 147
    assert np.array_equal(point[:147], np.array(baseline['parameters']))
    assert np.all(point[147:] == 0)
    residual, raw_jacobian, raw_profiles = raw.evaluate(point)
    checks = []

    def check(name, condition, value=None):
        checks.append(dict(name=name, passed=bool(condition), value=value))
        assert condition, (name, value)

    check('baseline objective replay', abs(residual @ residual - baseline['chi2']) < 1e-6,
          float(residual @ residual))
    models = [
        ('full_gls', GLSModel(free_shifts=True)),
        ('core01_gls', JointModel(free_shifts=True, threshold=.1, padding_km_s=0)),
        ('core03_gls', JointModel(free_shifts=True, threshold=.3, padding_km_s=1)),
    ]
    pattern = np.array([.1 if key in (2382, 2600) else 0. for key in raw.shift_keys])
    direction = pattern / .1
    cases = []
    saved = dict(common_parameters=point, illustrative_pattern_km_s=pattern)
    all_fisher = {}
    for name, model in models:
        rw, jacobian, profiles = model.evaluate(point)
        check(name + ': parameter labels common', model.labels == raw.labels)
        check(name + ': covariance kernel common',
              np.array_equal(model.kernel, models[0][1].kernel))
        check(name + ': forward profiles unchanged',
              all(np.array_equal(a, b) for a, b in zip(profiles, raw_profiles)))
        choleskies = model.cholesky if name == 'full_gls' else model.chols
        per_line = []
        # Independently check J^T C^-1 J for each line using a direct symmetric
        # positive-definite solve, without using the model's Cholesky whitener.
        for line, raw_line, sl, raw_sl, chol in zip(
                model.lines, raw.lines, model.slices, raw.slices, choleskies):
            native = line['idx'][line['good']]
            original = raw_line['idx'][raw_line['good']]
            take = np.searchsorted(original, native)
            check(f'{name} {line["key"]}: subset preserves native indices',
                  np.array_equal(original[take], native))
            direct_j = raw_jacobian[raw_sl][take]
            lag = np.abs(native[:, None] - native[None, :])
            covariance = np.zeros(lag.shape)
            inside = lag < len(model.kernel)
            covariance[inside] = model.kernel[lag[inside]]
            eigmin = float(np.linalg.eigvalsh(covariance)[0])
            check(f'{name} {line["key"]}: covariance positive definite', eigmin > 0, eigmin)
            check(f'{name} {line["key"]}: covariance reconstruction',
                  np.max(np.abs(chol @ chol.T - covariance)) < 2e-14)
            # Test all five shift columns plus two independent gas/continuum
            # combinations; full normal equations are intentionally avoided.
            rng = np.random.default_rng(260921 + line['key'])
            probe = np.column_stack([np.eye(model.npar)[:, 147:],
                                     rng.normal(size=(model.npar, 2))])
            x = direct_j @ probe
            direct_gram = x.T @ solve(covariance, x, assume_a='pos')
            whitened = jacobian[sl] @ probe
            rel = float(np.linalg.norm(direct_gram - whitened.T @ whitened) /
                        max(np.linalg.norm(direct_gram), 1.))
            check(f'{name} {line["key"]}: independent covariance solve', rel < 1e-12, rel)
            per_line.append(dict(line=line['key'], retained_pixels=int(len(native)),
                                 covariance_min_eigenvalue=eigmin))
        nuisance = jacobian[:, :147]
        shifts = jacobian[:, 147:]
        u, singular, vh = np.linalg.svd(nuisance, full_matrices=False)
        case_results = []
        for cutoff in CUTOFFS:
            rank = int(np.sum(singular > singular[0] * cutoff))
            basis = u[:, :rank]
            effective = shifts - basis @ (basis.T @ shifts)
            us, ss, vst = np.linalg.svd(effective, full_matrices=False)
            shift_rank = int(np.sum(ss > ss[0] * cutoff))
            fisher = effective.T @ effective
            covariance = (vst.T / ss**2) @ vst
            sigma = np.sqrt(np.diag(covariance))
            overlap = float(np.max(np.abs(basis.T @ effective)) /
                            max(np.max(np.abs(effective)), 1.))
            # QR is applied to the SAME retained physical nuisance subspace,
            # constructed from A V_r. This tests projection via a separate
            # decomposition without changing the singular cutoff convention.
            retained_design = nuisance @ vh[:rank].T
            q, _ = qr(retained_design, mode='economic')
            effective_qr = shifts - q @ (q.T @ shifts)
            qr_fisher = effective_qr.T @ effective_qr
            qr_relative = float(np.linalg.norm(qr_fisher - fisher) / np.linalg.norm(fisher))
            eig = np.linalg.eigvalsh(fisher)
            inverse_error = float(np.max(np.abs(fisher @ covariance - np.eye(5))))
            prefix = f'{name} cutoff={cutoff:g}'
            check(prefix + ': shift rank five', shift_rank == 5, shift_rank)
            check(prefix + ': projection orthogonality', overlap < 1e-11, overlap)
            check(prefix + ': effective information positive definite', eig[0] > 0, float(eig[0]))
            check(prefix + ': covariance inverse', inverse_error < 1e-10, inverse_error)
            check(prefix + ': independent QR projection', qr_relative < 2e-6, qr_relative)
            coherent_info = float(direction @ fisher @ direction)
            expected_local_noncentrality = float(pattern @ fisher @ pattern)
            info = dict(relative_svd_cutoff=cutoff, nuisance_rank=rank,
                        nuisance_singular_values=singular.tolist(),
                        effective_shift_rank=shift_rank, shift_singular_values=ss.tolist(),
                        fisher_per_km2_s2=fisher.tolist(), covariance_km2_s2=covariance.tolist(),
                        marginal_sigma_m_s=(1000*sigma).tolist(),
                        coherent_amplitude_information_per_km2_s2=coherent_info,
                        coherent_amplitude_sigma_m_s=1000/np.sqrt(coherent_info),
                        illustrative_100m_s_expected_local_noncentrality=expected_local_noncentrality,
                        projection_overlap=overlap, qr_fisher_relative_discrepancy=qr_relative)
            case_results.append(info)
            all_fisher[name, cutoff] = fisher
            saved[f'{name}_cutoff_{cutoff:g}_effective_shift_jacobian'] = effective
            saved[f'{name}_cutoff_{cutoff:g}_fisher'] = fisher
            saved[f'{name}_cutoff_{cutoff:g}_covariance'] = covariance
        cases.append(dict(case=name, retained_pixels=model.ndata, per_line=per_line,
                          mask_audit=getattr(model, 'mask_audit', []), cutoffs=case_results))
    for case in cases:
        for result, full in zip(case['cutoffs'], cases[0]['cutoffs']):
            cutoff = result['relative_svd_cutoff']
            result['coherent_information_retention_vs_full'] = (
                result['coherent_amplitude_information_per_km2_s2'] /
                full['coherent_amplitude_information_per_km2_s2'])
            result['marginal_variance_information_retention_vs_full'] = (
                (np.array(full['marginal_sigma_m_s']) /
                 np.array(result['marginal_sigma_m_s']))**2).tolist()
            loss_eigen = np.linalg.eigvalsh(all_fisher['full_gls', cutoff] -
                                          all_fisher[case['case'], cutoff])
            result['full_minus_mask_information_eigenvalues'] = loss_eigen.tolist()
            check(f'{case["case"]} {cutoff:g}: removing data does not increase information',
                  loss_eigen[0] > -1e-6, float(loss_eigen[0]))
    dependencies = [BASELINE, ROOT/'results/espresso_null/null_cross.npz',
                    ROOT/'results/noise_covariance/controls.json', Path(__file__),
                    Path(__file__).with_name('espresso_conventional_null.py'),
                    Path(__file__).with_name('espresso_empirical_gls.py'),
                    Path(__file__).with_name('espresso_saturation_controls.py'),
                    Path(__file__).with_name('joint_controls_fit.py')]
    result = dict(
        baseline_file=str(BASELINE.relative_to(ROOT)),
        baseline_chi2=baseline['chi2'], baseline_active_bounds=baseline['active_bounds'],
        common_parameters=point.tolist(), nuisance_parameter_count=147,
        shift_keys=raw.shift_keys, anchor=2374, primary_relative_svd_cutoff=1e-10,
        covariance_kernel=models[0][1].kernel.tolist(), continuum_variance_transferred=False,
        pattern_definition='Predeclared illustrative tangent signal: shifts 2382 and 2600 = +0.1 km/s relative to 2374, other three shifts = 0; gas/continuum tangent nuisance is projected out.',
        illustrative_pattern_km_s=pattern.tolist(), cases=cases,
        limitations=[
            'No nonlinear optimization, observed-residual score, discovery significance, or cosmic-time-law fit is performed.',
            'All designs use the exact same earlier null parameter vector, which is not the GLS optimum or masked-data optimum.',
            'The 147-dimensional nuisance tangent permits unconstrained infinitesimal displacements; 12 baseline parameters are at bounds, so these are not bound-aware confidence intervals.',
            'Weak gas directions depend on the SVD cutoff; all three cutoffs are reported in native parameter units without column rescaling.',
            'The coherent two-line pattern is illustrative and motivated by the earlier pattern, not an independent physical prediction or blind confirmatory signal.',
            'Information loss can accompany a reduced observed fit improvement, but cannot identify saturation as the cause or quantitatively explain a nonlinear likelihood-ratio change.',
            'The continuum-derived covariance shape is assumed transferable into absorption and does not include cross-line calibration or model-error covariance.',
        ],
        source_sha256={str(p.relative_to(ROOT)): digest(p) for p in dependencies},
        validation=dict(check_count=len(checks), all_passed=all(x['passed'] for x in checks), checks=checks))
    (OUT/'information.json').write_text(json.dumps(result, indent=2)+'\n')
    np.savez_compressed(OUT/'information_arrays.npz', **saved)
    central = [next(x for x in case['cutoffs'] if x['relative_svd_cutoff']==1e-10)
               for case in cases]
    lines = [
        '# 强线核心遮罩与局部速度信息：同一参数点的对照', '',
        '本诊断回答一个限定问题：遮掉强线核心后，约 100 m/s 的两线共同相对位移本身还剩多少可测信息？它不重新拟合，不使用观测残差来计算发现显著性，也不拟合宇宙时间规律。', '',
        '三个方案均使用此前 `null_cross.json` 的完全相同的 147 个气体与连续谱参数，所有额外相对位移设为零。两种遮罩固定来自同一个此前零位移模型。协方差沿用连续谱实测 ACF 的 Bartlett 截断核，且按真实像素索引重新构建，未把连续谱方差缩放迁移到吸收线。', '',
        '宽遮罩请求在阈值区域两侧各扩展 1 km/s；按照既有遮罩实现，在约 0.4 km/s 的格点上向外取整为每侧 3 个像素（约 1.2 km/s）。', '',
        '事先指定的示例信号为 Fe II 2382 与 2600 同时相对 2374 移动 +100 m/s，其余三条线相对位移为零。这个模式是说明性诊断，来自此前看到的强线模式，不是独立物理预测。信息是在允许气体与连续谱作局部线性调整后计算的。', '',
        '| 方案 | 像素数 | 100 m/s 示例的预期局部信号平方 | 信息保留率 | 同向位移幅度的局部误差 (m/s) |',
        '|---|---:|---:|---:|---:|',
    ]
    for case, stat in zip(cases, central):
        label={'full_gls':'完整 GLS','core01_gls':'F < 0.1 核心遮罩 + GLS','core03_gls':'F < 0.3 且扩展 1 km/s + GLS'}[case['case']]
        lines.append(f'| {label} | {case["retained_pixels"]} | {stat["illustrative_100m_s_expected_local_noncentrality"]:.4f} | {100*stat["coherent_information_retention_vs_full"]:.2f}% | {stat["coherent_amplitude_sigma_m_s"]:.3f} |')
    lines += ['', '上述信号平方是固定线性化与给定协方差下的预期非中心参数，不是实测 Δχ²，也不是显著性。幅度误差要求其余位移严格遵循所指定模式；下表则允许五个相对位移独立变化并同时边缘化，因此两类误差含义不同。', '',
              '| 方案 | 2260 | 2344 | 2382 | 2586 | 2600 |',
              '|---|---:|---:|---:|---:|---:|']
    for case, stat in zip(cases, central):
        lines.append('| '+case['case']+' | '+' | '.join(f'{x:.3f}' for x in stat['marginal_sigma_m_s'])+' |')
    lines += ['', '表中是五个位移各自的局部边缘化标准差，单位 m/s。完整的 5×5 信息矩阵、协方差和 10⁻⁸、10⁻¹⁰、10⁻¹² 奇异值阈值结果见 `results/mask_information/information.json`。', '',
              '这个对照允许说：更宽的核心遮罩同时移除了相当一部分位移信息，因而观测拟合改善下降不能自动被解读为已查明饱和造成了偏移。它不能证明全部下降都由信息损失解释；实际非线性拟合的最优参数会移动，未建模误差也可能改变。', '',
              f'数值检查：{len(checks)} 项全部通过，包括零位移目标函数重现、模型谱一致、真实索引上的协方差正定性、直接协方差求解与白化结果对照、SVD/QR 投影对照、投影正交性、信息矩阵正定与逆关系、删数据不增加有效信息。', '',
              '另有独立复算 `results/mask_information/independent_check.json`：使用不同 LAPACK SVD 驱动与 5×5 正定求解，三个方案的信号平方与边缘误差均与本结果一致；该复算还直接重建了阈值、像素扩展与原始有效掩罩，并检查源文件及共同参数未被修改。', '',
              '限制：该参数点有 12 个参数接近拟合边界；无约束的局部切空间投影不等于受边界约束的非线性置信区间。147 个气体与连续谱方向存在退化，阈值敏感性保留在输出中。未加入仪器标定、跨谱线系统误差或气体模型误差协方差。本结果仅为局部灵敏度说明。', '',
              '复算：`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python code/mask_information.py`（从本实验目录运行）。', '']
    (ROOT/'reports/mask_information_cn.md').write_text('\n'.join(lines))
    print(json.dumps(dict(check_count=len(checks), results=[dict(
        case=c['case'], ndata=c['retained_pixels'],
        expected_noncentrality=s['illustrative_100m_s_expected_local_noncentrality'],
        retention=s['coherent_information_retention_vs_full'],
        marginal_sigma_m_s=s['marginal_sigma_m_s']) for c,s in zip(cases,central)]), indent=2))


if __name__ == '__main__':
    main()
