#!/usr/bin/env python3
"""Descriptive integrals of actual SQUAD flux; no synthetic line or cosmic fit.

Run from any directory. Inputs are the previously archived, decoded SQUAD DR1
spectra and the archived NIST ASD Ritz wavelengths. NumPy and Matplotlib suffice.
The continuum, transition identities and integration regions are held fixed.
Flux-deficit centroids describe entire blended profiles, not transition shifts.
"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
C_KM_S = 299792.458
LABELS = ["Fe II 1608", "Fe II 2344", "Fe II 2374", "Fe II 2382", "Fe II 2586", "Fe II 2600"]
TARGETS = {
    "J051707-441055": {"alias": "HE 0515−4414", "z": 1.1508, "bounds": (-25., 110.)},
    "J034943-381030": {"alias": "Q0347−3819", "z": 3.025, "bounds": (-80., 35.)},
}
RNG_SEED = 20260920
MC_DRAWS = 5000


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def slice_profile(data, rest, z, bounds):
    """Piecewise-constant pixels, with exact partial widths at range edges.

    Optical velocity is a coordinate convention at the nominal absorber redshift.
    We use arithmetic pixel edges in wavelength, then its linear velocity map.
    No interpolation or clipping of absorption depths is performed.
    """
    wave = data["wavelength_vacuum_heliocentric_AA"]
    centre = rest * (1. + z)
    v = C_KM_S * (wave / centre - 1.)
    edge = np.empty(len(v) + 1)
    edge[1:-1] = (v[:-1] + v[1:]) / 2.
    edge[0] = v[0] - (v[1] - v[0]) / 2.
    edge[-1] = v[-1] + (v[-1] - v[-2]) / 2.
    left, right = np.maximum(edge[:-1], bounds[0]), np.minimum(edge[1:], bounds[1])
    valid = data["valid"] & (right > left)
    width = (right - left)[valid]
    # The first moment of a piecewise-constant partial pixel is at segment centre.
    segment_v = ((left + right) / 2.)[valid]
    return {
        "v": segment_v, "dv": width, "flux": data["flux_normalized"][valid],
        "error": data["error_normalized"][valid],
        "expected_fluctuation": data["expected_fluctuation_normalized"][valid],
        "pixel_index": np.flatnonzero(valid), "rest": rest,
        "coverage_fraction": float(width.sum() / (bounds[1] - bounds[0])),
    }


def integrals(p, continuum_scale=1.):
    """Signed absorption integral and first moment, without positive clipping."""
    f = p["flux"] / continuum_scale
    width, v = p["dv"], p["v"]
    deficit = (1. - f) * width
    area = deficit.sum()
    ew = area * p["rest"] / C_KM_S
    centroid = np.sum(v * deficit) / area
    g_ew = -width * p["rest"] / C_KM_S / continuum_scale
    g_centroid = -width * (v - centroid) / area / continuum_scale
    return float(ew), float(centroid), g_ew, g_centroid


def covariance_error(gradient, error, pixel_index, rho):
    # This is a hypothetical AR(1) covariance, NOT covariance measured from data.
    # Using original pixel indices preserves gaps in the native sampling.
    separation = np.abs(pixel_index[:, None] - pixel_index[None, :])
    correlation = rho ** separation
    weighted = gradient * error
    return float(np.sqrt(weighted @ correlation @ weighted))


def profile_class(target, label):
    if target == "J051707-441055":
        return "descriptive_profile_only"
    if label in ("Fe II 2586", "Fe II 2600"):
        return "excluded_low_snr"
    if label in ("Fe II 2344", "Fe II 2374", "Fe II 2382"):
        return "screen_only_telluric_risk"
    return "descriptive_profile_only"


def measure(data, target, info, atomic, rng):
    rest = float(atomic["ritz_vacuum_wavelength_A"])
    label = atomic["label"]
    lo, hi = info["bounds"]
    p = slice_profile(data, rest, info["z"], (lo, hi))
    ew, vc, ge, gv = integrals(p)
    se, sv = np.linalg.norm(ge * p["error"]), np.linalg.norm(gv * p["error"])
    good_expected = np.all(np.isfinite(p["expected_fluctuation"]) & (p["expected_fluctuation"] > 0))
    fe = p["expected_fluctuation"] if good_expected else p["error"]
    status = profile_class(target, label)
    centroid_reportable = bool(status != "excluded_low_snr" and abs(ew / se) >= 10)
    continuum_ew, continuum_v = [], []
    for scale in (0.99, 1.01):
        e, v, _, _ = integrals(p, scale)
        continuum_ew.append(e - ew)
        continuum_v.append(v - vc)
    range_ew, range_v = [], []
    variants = [(lo - 10, hi + 10), (lo + 10, hi - 10),
                (lo - 10, hi), (lo, hi + 10)]
    for bounds in variants:
        e, v, _, _ = integrals(slice_profile(data, rest, info["z"], bounds))
        range_ew.append(e - ew)
        range_v.append(v - vc)
    # Parametric propagation about measured flux; it validates only the quoted
    # diagonal noise propagation, not the physical source model or contamination.
    noise = rng.normal(size=(MC_DRAWS, len(p["flux"]))) * p["error"]
    deficit = (1. - p["flux"] - noise) * p["dv"]
    mc_ew = deficit.sum(axis=1) * rest / C_KM_S
    mc_v = (deficit @ p["v"]) / deficit.sum(axis=1)
    result = {
        "target": target, "alias": info["alias"], "z_abs_reference": info["z"],
        "transition": label, "status": status,
        "centroid_reportable_as_profile_descriptor": centroid_reportable,
        "rest_wavelength_NIST_Ritz_A": rest,
        "rest_wavelength_uncertainty_equivalent_m_s": float(atomic["rest_wavelength_uncertainty_equivalent_m_s"]),
        "v_min_km_s": lo, "v_max_km_s": hi, "n_pixels": len(p["v"]),
        "coverage_fraction": p["coverage_fraction"],
        "median_continuum_to_error": float(np.median(1. / p["error"])),
        "flux_min": float(np.min(p["flux"])), "flux_percentile05": float(np.quantile(p["flux"], 0.05)),
        "fraction_pixels_flux_below_0p1": float(np.mean(p["flux"] < 0.1)),
        "fraction_pixels_flux_below_3sigma": float(np.mean(p["flux"] <= 3*p["error"])),
        "rest_EW_mA": 1000*ew, "rest_EW_error_diagonal_mA": 1000*float(se),
        "rest_EW_error_expected_fluctuation_diagonal_mA": 1000*float(np.linalg.norm(ge * fe)),
        "expected_fluctuation_valid_all_selected_pixels": bool(good_expected),
        "EW_diagonal_signal_to_noise": float(ew/se),
        "centroid_km_s": vc if centroid_reportable else None,
        "centroid_error_diagonal_km_s": float(sv) if centroid_reportable else None,
        "centroid_error_expected_fluctuation_diagonal_km_s": float(np.linalg.norm(gv * fe)) if centroid_reportable else None,
        "rest_EW_continuum_1pct_max_change_mA": 1000*float(max(abs(np.array(continuum_ew)))),
        "centroid_continuum_1pct_max_change_km_s": float(max(abs(np.array(continuum_v)))) if centroid_reportable else None,
        "rest_EW_range_10kms_max_change_mA": 1000*float(max(abs(np.array(range_ew)))),
        "centroid_range_10kms_max_change_km_s": float(max(abs(np.array(range_v)))) if centroid_reportable else None,
        "rest_EW_error_MC_mA": 1000*float(np.std(mc_ew, ddof=1)),
        "centroid_error_MC_km_s": float(np.std(mc_v, ddof=1)) if centroid_reportable else None,
        "mc_EW_to_analytic_error_ratio": float(np.std(mc_ew, ddof=1)/se),
        "mc_centroid_to_analytic_error_ratio": float(np.std(mc_v, ddof=1)/sv) if centroid_reportable else None,
    }
    for rho in (0.3, 0.6):
        tag = str(rho).replace(".", "p")
        result[f"rest_EW_error_assumed_AR1_rho{tag}_mA"] = 1000*covariance_error(ge, p["error"], p["pixel_index"], rho)
        result[f"centroid_error_assumed_AR1_rho{tag}_km_s"] = covariance_error(gv, p["error"], p["pixel_index"], rho) if centroid_reportable else None
    return result, p


def make_figure(rows, profiles):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42})
    fig, axes = plt.subplots(3, 2, figsize=(12.3, 12.0))
    colours = ["#174F79", "#B45618", "#387A49", "#9D346D", "#826226", "#6D54A3"]
    for col, (target, info) in enumerate(TARGETS.items()):
        rr = [row for row in rows if row["target"] == target]
        ax = axes[0, col]
        displayed_flux = []
        for i, row in enumerate(rr):
            if row["status"] == "excluded_low_snr":
                continue
            p = profiles[(target, row["transition"])]
            displayed_flux.extend(p["flux"])
            linestyle = "--" if row["status"] == "screen_only_telluric_risk" else "-"
            ax.step(p["v"], p["flux"], where="mid", color=colours[i], lw=.9,
                    ls=linestyle, label=row["transition"])
        ax.axhline(1., color="gray", lw=.7, ls=":")
        ax.set(xlabel="Optical velocity at nominal absorber redshift (km/s)",
               ylabel="Actual normalized flux", ylim=(min(-.1, min(displayed_flux)-.03), max(1.2, max(displayed_flux)+.03)),
               title=f"{info['alias']}  ·  z_abs = {info['z']}")
        ax.legend(fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(.5, -.20), frameon=False)
        ax.grid(alpha=.15)
        x = np.arange(6)
        ax = axes[1, col]
        for i, row in enumerate(rr):
            marker = "x" if row["status"] == "excluded_low_snr" else ("s" if row["status"] == "screen_only_telluric_risk" else "o")
            ax.errorbar(i, row["rest_EW_mA"], yerr=row["rest_EW_error_diagonal_mA"],
                        color=colours[i], marker=marker, capsize=3, ms=6, lw=1.1)
            ax.text(i, row["rest_EW_mA"]+row["rest_EW_error_diagonal_mA"]+10,
                    f"{row['rest_EW_mA']:.1f}", ha="center", fontsize=8)
        ax.set(xticks=x, xticklabels=[s.split()[-1] for s in LABELS],
               xlabel="Fe II nominal rest wavelength (Å)", ylabel="Signed rest equivalent width (mÅ)")
        ax.set_title("Measured finite-window absorption integrals", fontsize=10)
        ax.grid(axis="y", alpha=.18)
        ax = axes[2, col]
        for i, row in enumerate(rr):
            if row["centroid_km_s"] is None:
                continue
            y = row["centroid_km_s"]
            ax.errorbar(i, y, yerr=row["centroid_range_10kms_max_change_km_s"],
                        color="#AAB0BB", lw=6, alpha=.7, capsize=0)
            ax.errorbar(i, y, yerr=row["centroid_continuum_1pct_max_change_km_s"],
                        color="#D2A35B", lw=3, alpha=.8, capsize=0)
            marker = "s" if row["status"] == "screen_only_telluric_risk" else "o"
            ax.errorbar(i, y, yerr=row["centroid_error_diagonal_km_s"],
                        color=colours[i], marker=marker, capsize=4, ms=5, lw=1.3)
        ax.set(xticks=x, xticklabels=[s.split()[-1] for s in LABELS],
               xlabel="Fe II nominal rest wavelength (Å)",
               ylabel="Flux-deficit centroid (km/s)")
        ax.set_title("Profile descriptor; not a transition-frequency shift", fontsize=10)
        ax.grid(axis="y", alpha=.18)
        if col == 1:
            ax.text(.99, .04, "Squares / dashed: telluric risk\n2586, 2600 centroids excluded (low S/N)",
                    ha="right", va="bottom", transform=ax.transAxes, fontsize=8, color="#864216")
    fig.suptitle("Existing-data analysis: absorption strength, shape and robustness", fontsize=16, y=.985)
    fig.text(.5, .948, "Real SQUAD DR1 flux · NIST Ritz wavelengths · fixed exploratory integration ranges · no cosmic fit",
             ha="center", fontsize=10, color="#435366")
    fig.text(.5, .055, "Bottom: coloured bars = conditional diagonal 1σ; ochre = ±1% continuum sensitivity; grey = ±10 km/s range sensitivity.", ha="center", fontsize=8.5)
    fig.text(.5, .035, "High-z red windows are screening diagnostics. Native-pixel covariance, tellurics and detailed gas components remain unmodelled.", ha="center", fontsize=8.5)
    fig.subplots_adjust(left=.08, right=.98, bottom=.115, top=.903, hspace=.63, wspace=.26)
    for suffix in ("png", "pdf"):
        fig.savefig(ROOT / f"figures/observed_profile_diagnostics.{suffix}", dpi=180, facecolor="white")
    plt.close(fig)


def write_report(rows, summary):
    table = ["| 目标 | 谱线 | 静止系等效宽度 mÅ | 吸收质心 km/s | 连续谱 1% 最大改变量 km/s | 窗口改动最大改变量 km/s | 使用状态 |",
             "|---|---:|---:|---:|---:|---:|---|"]
    for r in rows:
        centroid = "不报告" if r["centroid_km_s"] is None else f"{r['centroid_km_s']:.3f} ± {r['centroid_error_diagonal_km_s']:.3f}"
        con = "—" if r["centroid_km_s"] is None else f"{r['centroid_continuum_1pct_max_change_km_s']:.3f}"
        ran = "—" if r["centroid_km_s"] is None else f"{r['centroid_range_10kms_max_change_km_s']:.3f}"
        state = {"descriptive_profile_only": "描述性测量", "screen_only_telluric_risk": "水汽风险，仅筛查", "excluded_low_snr": "低信噪，排除质心"}[r["status"]]
        table.append(f"| {r['alias']} | {r['transition'].split()[-1]} | {r['rest_EW_mA']:.2f} ± {r['rest_EW_error_diagonal_mA']:.2f} | {centroid} | {con} | {ran} | {state} |")
    text = """# 两份真实档案光谱的描述性分析

本轮直接使用已下载 UVES SQUAD DR1 合并光谱中的真实归一化通量，测量有限窗口的等效宽度、吸收缺失加权质心和深吸收比例；没有注入谱线，也没有拟合宇宙时间趋势。它把“公开数据可以用”推进为实际可复算的测量与质量诊断。

## 定义与范围

采用已归档 NIST ASD 的真空 Ritz 波长重新计算光学速度坐标 `v = c[λ/(λ₀(1+z_abs))−1]`。日心真空波长与原生像素保持不变。积分边界穿过像素时按边界切分像素宽度；没有重采样，没有把负吸收（通量高于连续谱）裁成零。静止系等效宽度为 `W_r = (λ₀/c) Σ(1−F_i)Δv_i`；质心为 `v_bar = Σv_i(1−F_i)Δv_i / Σ(1−F_i)Δv_i`。边界像素的一阶矩使用截取片段中点。

HE 0515−4414 使用 `z_abs=1.1508`、`v∈[−25,110] km/s`；Q0347−3819 使用 `z_abs=3.025`、`v∈[−80,35] km/s`。这些是查看真实谱形后选定、在各跃迁间共用的探索性窗口，不是事前注册的盲选择，也不是整个吸收系统的总等效宽度。窗口选择旨在截取主要吸收复合体，因此不能将结果与采用整个系统积分的文献值直接相比。

## 实测结果

表中的 ± 仅为原文件归一化误差数组、固定连续谱和固定波长尺度、忽略像素协方差时的条件统计 1σ；不包含污染、分量结构、连续谱或波长校准系统误差。所有谱窗覆盖分数均在机器可读结果中检查。

TABLE

## 误差与稳健性

1. 对等效宽度做线性误差传播，对比值型质心做一阶梯度传播。用 5,000 次围绕真实测量通量的高斯噪声重抽样检查传播计算；这不是重新生成一条模型谱线，也不是物理模型验证。
2. 同时使用源文件提供的 `expected_fluctuation` 数组重复对角传播。SQUAD 文档说明该数组考虑重分散对像素间预期波动的影响，在谱线拟合中有用；它仍不能代替完整协方差矩阵。
3. 以假定 `Cov(F_i,F_j)=σ_iσ_jρ^|i−j|`、`ρ=0.3,0.6` 展示相关像素可能带来的误差变化。这两个相关系数是敏感性情景，不是从这些谱窗估计的实测协方差。
4. 将连续谱分别乘以 `0.99,1.01`、使归一化通量相应除以该因子，记录相对基准的最大绝对改变量。这是人为设置的 1% 扰动，不是已测得的连续谱不确定度，不加入统计 1σ。
5. 用扩展两端、收缩两端、只扩左端、只扩右端四种 ±10 km/s 范围变化评估窗口敏感性。报告的是最大绝对改变量，不是概率置信区间。范围变化会真实地增加或删除气体分量。
6. 深吸收指标为最小通量、通量第 5 百分位、`F<0.1` 像素比例和 `F≤3σ` 比例；它们提示饱和风险，不能单独证明未解析饱和。没有从近零/负通量计算并宣称可靠的光学深度或柱密度。

## 对假说检验的实际意义

HE 六个谱窗的条件统计质心误差约 18–117 m/s，而 1% 连续谱扰动导致约 99–439 m/s 的变化，窗口敏感性达到约 0.81–1.66 km/s。后两项不是实测误差条，却说明仅有很高的每像素信噪比还不足以支持极小频移的解释。HE 2382 和 2600 Å 谱窗分别有约 36.2% 和 31.7% 的像素低于归一化通量 0.1，属于应优先检查饱和影响的深吸收轮廓。

同一离子的不同跃迁已经出现可见的形状与吸收质心差异。这些差异可由复杂气体分量、不同吸收强度、饱和、混叠和仪器响应产生；本轮尚未分解各因素的贡献。仅靠这组描述量不能将其识别为原子跃迁频移，更不能换算成精细结构常数变化。质心还取决于参考红移，HE 的约 34–36 km/s 与 Q0347 的约 −20 km/s 不是绝对宇宙频移残差。低红移样例的高信噪适合开发分量联合拟合与系统误差检验；高红移样例说明大气吸收和波长覆盖是实际限制。

表中的统计误差没有加入 NIST Ritz 波长的不确定度（所用 Fe II 行约 13–15 m/s），其相关性以及星光波长标定的系统误差也未在公开合并谱中确定，因此没有构成 ppm 级跃迁比较的完整误差预算。

Q0347 的 2344/2374/2382 Å 对应观测波长约 9435–9591 Å，具有水汽污染风险，此处数字仅用于筛查，不能视为已经校正的大气外 Fe II 测量。2586/2600 Å 附近低信噪只保留带符号的窗口积分和质量指标，不给科学可用的质心；大等效宽度也不能使不可靠谱窗自动变得可用。

这两条视线是便利先导样本，不能用于声称红移演化、`1/ln²(t)` 优于其他形式或物理黎曼零点截断。数据分析所支持的是：现在可以开展实际吸收谱分析，但未来检验需固定原子物理响应、独立筛查与校准，随后再用更大且同质的样本检验时间规律。

## 复现与来源

运行 `python code/measure_observed_profiles.py`。结果为 `results/observed_profile_metrics.csv`、同名 JSON 和 `figures/observed_profile_diagnostics.pdf/.png`。JSON 保留输入 SHA-256、方法定义、噪声模型和 12 个谱窗完整诊断。脚本只依赖已经归档的输入，运行时不联网。

- [ESO SQUAD DR1 发布](https://www.eso.org/sci/publications/announcements/sciann17192.html)；本地原作者格式说明 `data/raw/Notes_FITS_Files.txt` 解释误差与 expected-fluctuation 数组。
- [NIST ASD 谱线说明](https://www.physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html)；具体 Ritz 值与引用代码见 `data/atomic/selected_transitions.csv`。
- [Savage & Sembach (1991), DOI 10.1086/170498](https://doi.org/10.1086/170498) 讨论吸收谱形、表观光学深度及其解释；本轮仅做通量积分，未进行未解析饱和修正或柱密度估计。
""".replace("TABLE", "\n".join(table))
    (ROOT / "reports/observed_profiles_cn.md").write_text(text, encoding="utf-8")


def main():
    for name in ("results", "figures", "reports"):
        (ROOT / name).mkdir(exist_ok=True)
    atomic_path = ROOT / "data/atomic/selected_transitions.csv"
    atomic = {r["label"]: r for r in csv.DictReader(atomic_path.open())}
    rng = np.random.default_rng(RNG_SEED)
    rows, profiles = [], {}
    source_hashes = {str(atomic_path.relative_to(ROOT)): sha256(atomic_path)}
    for target, info in TARGETS.items():
        path = ROOT / f"data/processed/{target}_squad_dr1.npz"
        source_hashes[str(path.relative_to(ROOT))] = sha256(path)
        with np.load(path, allow_pickle=False) as f:
            data = {key: f[key] for key in f.files}
        for label in LABELS:
            row, profile = measure(data, target, info, atomic[label], rng)
            rows.append(row)
            profiles[(target, label)] = profile
    with (ROOT / "results/observed_profile_metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "analysis_date": "2026-09-20", "data_kind": "actual observed normalized SQUAD DR1 flux",
        "n_spectra": 2, "n_windows": 12, "n_mc_draws": MC_DRAWS, "rng_seed": RNG_SEED,
        "input_sha256": source_hashes,
        "selection": "Exploratory visually selected shared velocity windows; two convenience targets, not population sample",
        "velocity_convention": "optical c*(lambda/(NIST_Ritz_rest*(1+z_abs))-1); vacuum heliocentric wavelength",
        "integration": "piecewise-constant flux, midpoint wavelength-cell edges, exact partial widths and first moments at bounds; signed deficits, no clipping",
        "uncertainty": "conditional diagonal analytic 1sigma; independent Gaussian Monte Carlo only checks propagation",
        "expected_fluctuation": "alternative diagonal input, not full covariance; source documentation recommends it for model fitting",
        "covariance_sensitivity": "AR1 rho=.3,.6 assumed scenarios, not estimated covariance",
        "continuum_sensitivity": "input continuum multiplied by .99 or 1.01; max absolute change, not 1sigma",
        "window_sensitivity": "lo-10/hi+10, lo+10/hi-10, lo-10/hi, lo/hi+10 in km/s; max absolute change",
        "centroid_exclusion": "last two high-z windows excluded for low SNR; red high-z lines remain screen-only because tellurics uncorrected",
        "inference_limit": "profile centroids are not intrinsic transition shifts; no alpha, time-law, or zero-cutoff inference",
        "all_window_coverage_complete": all(abs(r["coverage_fraction"]-1)<1e-10 for r in rows),
        "rows": rows,
    }
    (ROOT / "results/observed_profile_metrics.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False)+"\n")
    make_figure(rows, profiles)
    write_report(rows, summary)
    for row in rows:
        print(row["target"], row["transition"], "EW", round(row["rest_EW_mA"], 3), "mA", "centroid", row["centroid_km_s"], "km/s", row["status"])


if __name__ == "__main__":
    main()
