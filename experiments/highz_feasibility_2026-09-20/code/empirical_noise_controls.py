#!/usr/bin/env python3
"""Nearby continuum-control ACF audit, followed by fixed local covariance tests.

This is exploratory: the spectrum and previous absorption residuals were already
inspected. No interval is selected by its ACF or by whether it reduces a shift.
All fixed intervals, exclusions, threshold variants and pixel ownership are saved.
No per-pixel residual clipping is performed; covariance transfer into absorption
pixels is a conditional sensitivity, not a recovered instrumental covariance.
"""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import json, hashlib, csv
import numpy as np
from scipy.linalg import toeplitz, cholesky, solve_triangular
from espresso_conventional_null import ROOT, RAW, Model, load_sources, C, ZREF
from espresso_null_local_controls import project_controls

OUT=ROOT/'results/noise_covariance'
PROC=ROOT/'data/processed'
REPORT=ROOT/'reports/noise_covariance_controls_cn.md'
LOWERS=[-1320,-1100,-880,180,400,620]
WIDTH=220
MAXLAG=10
SEED=260920375

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def detrend(idx, velocity, flux, error, degree=1):
    x=(velocity-velocity.mean())/(np.ptp(velocity) or 1)
    X=np.polynomial.polynomial.polyvander(x,degree)
    coeff=np.linalg.lstsq(X/error[:,None],flux/error,rcond=None)[0]
    trend=X@coeff
    return (flux-trend)/error, trend

def moments(idx,r):
    """Lag products only for consecutive retained source pixels, never gaps."""
    cross=np.zeros(MAXLAG+1);pairs=np.zeros(MAXLAG+1,dtype=int)
    cross[0]=r@r;pairs[0]=len(r)
    for lag in range(1,MAXLAG+1):
        good=idx[lag:]-idx[:-lag]==lag
        cross[lag]=np.sum(r[lag:][good]*r[:-lag][good]);pairs[lag]=sum(good)
    return cross,pairs

def aggregate(rows, bootstrap=0):
    if not rows:return None
    cross=np.array([x['cross'] for x in rows]);pairs=np.array([x['pairs'] for x in rows])
    covariance=cross.sum(0)/pairs.sum(0)
    result=dict(interval_count=len(rows),pixel_count=int(pairs[:,0].sum()),
        variance=float(covariance[0]),rms=float(np.sqrt(covariance[0])),
        covariance=covariance.tolist(),acf=(covariance/covariance[0]).tolist(),
        lag_pair_counts=pairs.sum(0).tolist())
    if bootstrap:
        rng=np.random.default_rng(SEED)
        draws=rng.integers(len(rows),size=(bootstrap,len(rows)))
        cov=cross[draws].sum(1)/pairs[draws].sum(1)
        acf=cov/cov[:,0,None]
        result['interval_cluster_bootstrap']=dict(draws=bootstrap,seed=SEED,
            variance_95_percentile_interval=np.quantile(cov[:,0],[.025,.975]).tolist(),
            acf_95_percentile_intervals=np.quantile(acf,[.025,.975],axis=0).T.tolist(),
            caveat='Resamples selected disjoint wavelength intervals as clusters. It is a conditional sampling interval, not a confidence interval including selection or wavelength-dependent systematics.')
    return result

def screen(row,threshold=8.,rms_limit=1.25,min_fraction=.9,median_bounds=(.97,1.03)):
    reasons=[]
    if row['retained_fraction']<min_fraction: reasons.append(f'less_than_{100*min_fraction:g}_percent_unique_outside_author_fit_regions')
    if not (median_bounds[0]<=row.get('median_flux',0)<=median_bounds[1]):reasons.append(f'median_flux_outside_{median_bounds[0]:g}_to_{median_bounds[1]:g}')
    if row.get('max_abs_20pixel_bin_score',np.inf)>threshold:reasons.append('broad_residual_bin_threshold')
    if row.get('rms',np.inf)>rms_limit:reasons.append('standardized_residual_rms_threshold')
    return reasons

def controls():
    comp,atoms,regions,wave,a,dv=load_sources()
    # Exclude all published species/regions, not just the six lines in our fit.
    known=np.zeros(len(wave),bool);published=[]
    for path in sorted(RAW.glob('fit_[lcr]_iso.f13')):
        for line in path.read_text().splitlines():
            if line.startswith('hes'):
                tok=line.split();lo,hi=map(float,tok[2:4]);known|=(wave>=lo)&(wave<=hi)
                published.append(dict(file=path.name,lo_AA=lo,hi_AA=hi))
    assigned=np.zeros(len(wave),bool);rows=[];arrays={}
    for k,atom in atoms.items():
        ref=np.sum(atom[:,0]*atom[:,1])/atom[:,1].sum()
        v=C*np.log(wave/(ref*(1+ZREF)))
        for lo in LOWERS:
            candidate=np.flatnonzero((v>=lo)&(v<lo+WIDTH))
            valid=(a[4,candidate]==1)&np.isfinite(a[0,candidate])&np.isfinite(a[1,candidate])&(a[1,candidate]>0)
            unique=~assigned[candidate];assigned[candidate]=True
            good=valid&unique&~known[candidate]
            idx=candidate[good]
            name=f'{k}_{lo:+d}_{lo+WIDTH:+d}'
            row=dict(name=name,line=k,side='red' if lo>0 else 'blue',velocity_lo_km_s=lo,velocity_hi_km_s=lo+WIDTH,
                candidate_pixel_count=len(candidate),retained_pixel_count=len(idx),retained_fraction=len(idx)/len(candidate),
                author_fit_region_overlap=int(known[candidate].sum()),duplicate_pixels=int((~unique).sum()),invalid_pixels=int((~valid).sum()))
            if len(idx)>=50:
                r,trend=detrend(idx,v[idx],a[0,idx],a[1,idx])
                # Non-overlapping fixed 20-source-pixel bins. Keep only complete
                # bins, so missing pixels cannot create artificial adjacency.
                bins=[]
                for start in range(int(candidate[0]),int(candidate[-1])+1,20):
                    take=(idx>=start)&(idx<start+20)
                    if sum(take)==20:bins.append(float(np.sum(r[take])/np.sqrt(20)))
                cross,pairs=moments(idx,r)
                row.update(median_flux=float(np.median(a[0,idx])),rms=float(np.sqrt(np.mean(r*r))),
                    max_abs_20pixel_bin_score=max(map(abs,bins)) if bins else float('inf'),
                    complete_20pixel_bin_count=len(bins),cross=cross.tolist(),pairs=pairs.tolist(),
                    acf=(cross/pairs/(cross[0]/pairs[0])).tolist())
                # Alternate continuum degrees use identical pixels and screening.
                row['alternative_detrend_moments']={}
                for degree in [0,2]:
                    rd,_=detrend(idx,v[idx],a[0,idx],a[1,idx],degree)
                    cd,pd=moments(idx,rd)
                    row['alternative_detrend_moments'][str(degree)]={'cross':cd.tolist(),'pairs':pd.tolist()}
                arrays[name+'_source_index']=idx;arrays[name+'_standardized_residual']=r
                arrays[name+'_wavelength_AA']=wave[idx];arrays[name+'_flux']=a[0,idx]
                arrays[name+'_error']=a[1,idx];arrays[name+'_continuum_trend']=trend
            row['screen_exclusions']=screen(row);row['selected']=not row['screen_exclusions']
            rows.append(row)
    selected=[r for r in rows if r['selected']]
    all_adequate=[r for r in rows if r['retained_fraction']>=.9 and 'cross' in r]
    variants={}
    for threshold in [6.,8.,10.,12.]:
        for rmslim in [1.1,1.25,1.5]:
            chosen=[r for r in rows if not screen(r,threshold,rmslim)]
            variants[f'bin{threshold:g}_rms{rmslim:g}']=aggregate(chosen)
    result=dict(source_sha256=sha(RAW/'hes0515m4414.fits'),pixel_spacing_km_s=dv,
        selection='Exploratory fixed 220 km/s intervals at six specified offsets per line; all author fitting regions and duplicate source pixels excluded; no per-pixel residual clipping.',
        selection_thresholds=dict(minimum_retained_fraction=.9,median_flux=[.97,1.03],maximum_abs_20pixel_bin_score=8.,maximum_standardized_rms=1.25),
        bin_score_note='sum of standardized residuals / sqrt(20), used only as a deterministic broad-feature screen, not a Gaussian significance in correlated noise.',
        published_regions_excluded=published,intervals=rows,primary=aggregate(selected,bootstrap=10000),
        unfiltered_adequate_intervals=aggregate(all_adequate),threshold_sensitivities=variants,
        geometric_retention_sensitivities={str(frac):aggregate([r for r in rows if not screen(r,min_fraction=frac)]) for frac in [.5,.7,.8,.9,1.]},
        median_gate_sensitivities={f'{low}_{high}':aggregate([r for r in rows if not screen(r,median_bounds=(low,high))]) for low,high in [(.95,1.05),(.97,1.03),(.99,1.01)]},
        by_side={side:aggregate([r for r in selected if r['side']==side]) for side in ['red','blue']},
        by_transition={str(k):aggregate([r for r in selected if r['line']==k]) for k in atoms},
        detrending_sensitivities={str(degree):aggregate([r['alternative_detrend_moments'][str(degree)] for r in selected]) for degree in [0,2]},
        caveats=['Control spectra may contain weak undetected absorption, continuum errors, extraction artifacts and residual sky features.',
            'No covariance between different wavelength intervals is estimated.',
            'Flux-dependent errors and clipping can differ inside absorption; continuum scaling is not automatically transferable.',
            'The ACF does not identify the physical origin of correlation or validate an AR(1) process.'])
    np.savez_compressed(PROC/'espresso_noise_controls.npz',**arrays)
    (OUT/'controls.json').write_text(json.dumps(result,indent=2)+'\n')
    with (OUT/'intervals.csv').open('w') as f:
        fields=['name','line','side','velocity_lo_km_s','velocity_hi_km_s','candidate_pixel_count','retained_pixel_count','author_fit_region_overlap','duplicate_pixels','invalid_pixels','median_flux','rms','max_abs_20pixel_bin_score','selected','screen_exclusions']
        w=csv.DictWriter(f,fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
    return result

def whiten_segments(model,r,J,kind,rho,variance,acf):
    rw=np.empty_like(r);jw=np.empty_like(J);mineigen=1.
    for line,sl in zip(model.lines,model.slices):
        idx=line['idx'][line['good']];starts=np.r_[0,np.flatnonzero(np.diff(idx)>1)+1,len(idx)]
        for beg,end in zip(starts[:-1],starts[1:]):
            n=end-beg;lags=np.arange(n)
            if kind=='ar1': first=rho**lags
            elif kind=='ma1':
                first=np.zeros(n);first[0]=1
                if n>1:first[1]=rho
            elif kind=='tapered_acf10':
                first=np.zeros(n);m=min(n,len(acf));first[:m]=acf[:m]*(1-np.arange(m)/len(acf))
            else: first=np.r_[1.,np.zeros(n-1)]
            cov=toeplitz(first)*variance
            chol=cholesky(cov,lower=True)
            mineigen=min(mineigen,float(np.linalg.eigvalsh(cov)[0]))
            dest=slice(sl.start+beg,sl.start+end)
            rw[dest]=solve_triangular(chol,r[dest],lower=True)
            jw[dest]=solve_triangular(chol,J[dest],lower=True)
    return rw,jw,mineigen

def local_sensitivity(control):
    path=ROOT/'results/espresso_null/null_cross.json'
    old=json.loads(path.read_text());model=Model(free_shifts=True,**old['configuration'])
    p=model.initial()[0];lookup=dict(zip(old['labels'],old['parameters']))
    for i,key in enumerate(model.labels):p[i]=lookup.get(key,p[i])
    r,J,_=model.evaluate(p)
    assert abs(r@r-old['chi2'])<1e-6
    n_nuisance=model.shift_offset
    acf=np.array(control['primary']['acf']);rho=float(acf[1]);variance=control['primary']['variance']
    cases=[('diagonal_nominal','diagonal',0.,1.),('ar1_030_nominal','ar1',.3,1.),('ar1_040_nominal','ar1',.4,1.),('ar1_050_nominal','ar1',.5,1.)]
    for kind in ['ar1','ma1','tapered_acf10']:
        for scale,var in [('nominal',1.),('control_scaled',variance)]:cases.append((kind+'_'+scale,kind,rho,var))
    results=[]
    for name,kind,rhov,var in cases:
        rw,jw,mineig=whiten_segments(model,r,J,kind,rhov,var,acf)
        for cutoff in [1e-8,1e-10,1e-12]:
            summary,_=project_controls(jw,rw,n_nuisance,cutoff)
            # Remove conditional p-value to discourage an inferential interpretation.
            summary.pop('conditional_chi2_survival',None)
            summary.update(name=name,covariance_family=kind,rho1=rhov if kind!='diagonal' else 0.,variance_factor=var,minimum_covariance_eigenvalue=mineig,
                interpretation='Fixed-baseline unconstrained local nuisance projection under an assumed covariance. No nonlinear refit; no bounds, model selection, covariance-estimation uncertainty or discovery significance accounted for.')
            results.append(summary)
    out=dict(baseline=str(path.relative_to(ROOT)),baseline_sha256=sha(path),baseline_chi2=old['chi2'],
        baseline_nuisance_parameter_count=n_nuisance,shift_lines=model.shift_keys,reference_line=2374,
        covariance_note='Reset at every masked gap and transition. MA(1) uses only measured lag 1. Tapered ACF uses measured lags 0-10 times Bartlett weights 1-lag/11, then zero beyond lag 10. No covariance family is asserted to be physically correct.',
        empirical_rho1=rho,empirical_variance=variance,cases=results)
    (OUT/'local_shift_covariance_sensitivity.json').write_text(json.dumps(out,indent=2)+'\n')
    return out

def write_report(c,s):
    p=c['primary'];ci=p['interval_cluster_bootstrap']['acf_95_percentile_intervals'][1]
    lines=['# ESPRESSO 邻近连续谱控制与像素协方差敏感性','',
        '**这里实际测量了吸收区之外的像素相关；将该相关移入吸收模型的部分只是固定基线的局部敏感性分析，没有完成新的非线性重拟合。**','',
        '## 数据与选择','',
        '来源为已归档的 ESPRESSO 最终合并光谱。围绕每条 Fe II 参考波长，固定选取 6 个宽 220 km/s 的区间：[-1320,-1100]、[-1100,-880]、[-880,-660]、[180,400]、[400,620]、[620,840] km/s。处理按波长标签顺序进行，重叠原始像素只归入首个区间；排除原作者全部三个分区的正式拟合波长区间和 FITS 无效像素。',
        '全部 36 个候选区间的通过和排除原因保存在 intervals.csv。保留比例至少 90%、中位归一化通量在 0.97–1.03、线性连续谱去趋势后的标准化残差 RMS 不超过 1.25、固定 20 像素分箱的残差和除以 sqrt(20) 不超过 8。最后一项只是确定性的宽特征筛选，不能作为相关噪声下的 8σ 判据。没有逐像素残差裁剪，也没有按自相关或最终位移大小挑区间。',
        '这不是预注册或盲分析：此前已看过数据，阈值在初步检查后制定，因此同时给出宽松/严格阈值、红侧/蓝侧、逐跃迁和常数/二次连续谱敏感性。','',
        '## 实测结果','',
        f'- 通过筛选 {p["interval_count"]} 个互不重复区间，共 {p["pixel_count"]} 个像素。',
        f'- 残差相对于文件提供误差的 RMS = {p["rms"]:.5f}，方差因子 = {p["variance"]:.5f}。',
        f'- 邻像素相关 r(1) = {p["acf"][1]:.5f}；按整个区间重抽样的条件性 95% 百分位区间 [{ci[0]:.5f}, {ci[1]:.5f}]。',
        f'- r(2) = {p["acf"][2]:.5f}，r(3) = {p["acf"][3]:.5f}。全部 1–10 阶结果见 controls.json。',
        f'- 红侧 {c["by_side"]["red"]["interval_count"]} 区间的 r(1) = {c["by_side"]["red"]["acf"][1]:.5f}；蓝侧 {c["by_side"]["blue"]["interval_count"]} 区间的 r(1) = {c["by_side"]["blue"]["acf"][1]:.5f}。',
        '- 12 种分箱/RMS 阈值组合均保留同样 25 个区间：主要排除已由已知拟合区间或重复像素的几何规则完成。这并不证明区间位置、像素归属或其他选择没有影响；另列几何保留比例和中位通量门限的敏感性。',
        '- 方差采用连续谱去趋势后的残差均方，未另作自由度修正；常数/二次去趋势结果单列。对相关噪声，不能机械使用独立样本的 n−2 修正。',
        '- 此相关确实存在于邻近连续谱，不能全部归因于吸收气体模型残差。但其成因可能包括重采样、连续谱、弱吸收和提取过程；连续谱测量不能恢复原始曝光的完整协方差。','',
        '## 吸收模型的条件性局部复核','',
        '保留 null_cross 最优点和同一批 2,931 个像素，在白化后用直接 SVD 投影掉 147 个气体/连续谱方向，再计算 5 个相对位移方向的局部分数。表中数字不是新非线性最小 χ²，也不是发现显著性；界限约束、气体模型搜索和协方差估计不确定性未纳入。主表用相对奇异值截断 1e-10，完整输出另外保存 1e-8 和 1e-12。','',
        '| 固定协方差敏感性 | 方差因子 | 局部位移分数 | 2382 位移 (m/s) | 2600 位移 (m/s) |','|---|---:|---:|---:|---:|']
    for row in s['cases']:
        if row['relative_singular_value_cutoff']!=1e-10:continue
        i2382=s['shift_lines'].index(2382);i2600=s['shift_lines'].index(2600)
        lines.append(f'| {row["name"]} | {row["variance_factor"]:.4f} | {row["local_score_delta_chi2"]:.3f} | {row["local_shift_estimates_km_s"][i2382]*1000:.2f} | {row["local_shift_estimates_km_s"][i2600]*1000:.2f} |')
    lines+=['',
        'AR(1) 只强制 r(k)=r(1)^k；实测高阶 ACF 不必符合它。MA(1) 只保留邻像素相关；tapered_acf10 对实测 0–10 阶相关乘 Bartlett 权重 1−k/11，随后截断。所有协方差均逐连续有效像素段应用，并经 Cholesky 正定性检查。',
        '只加相关并保持原误差，通常会扩大平滑位移方向的不确定性；若同时把文件误差缩小到连续谱实测 RMS，这种下降可能被部分抵消。因此必须同时展示 nominal 和 control_scaled。后者尤其不能视作最终误差校正，因为深吸收区的光子计数、背景、裁剪和合并权重都可能不同。','',
        '## 文献边界','',
        '原始论文报告用 UVES_popler 合并 17 次曝光，进行了自动和人工剔除、手动连续谱处理，并用改变输出色散后重新拟合来评估重采样效应。本项控制是额外的最终合并谱检验，不能代替该曝光级重分析。[Murphy et al. 2022, §§2.3, 4.2.2](https://arxiv.org/abs/2112.05819)。','',
        '复现：`OPENBLAS_NUM_THREADS=1 python code/empirical_noise_controls.py`（本实验目录下）。全部选择、数字和协方差族保存在 results/noise_covariance；控制像素及趋势保存在 data/processed/espresso_noise_controls.npz。','']
    REPORT.write_text('\n'.join(lines))

def main():
    OUT.mkdir(parents=True,exist_ok=True);PROC.mkdir(parents=True,exist_ok=True)
    c=controls();print(json.dumps(c['primary'],indent=2),flush=True)
    s=local_sensitivity(c);write_report(c,s)
    dependencies=[Path(__file__).with_name('espresso_conventional_null.py'),Path(__file__).with_name('espresso_null_local_controls.py'),RAW/'MM_VPFIT_2013-11-10.dat',RAW/'fit_r_iso.f18',*sorted(RAW.glob('fit_[lcr]_iso.f13'))]
    (OUT/'provenance.json').write_text(json.dumps(dict(script_sha256=sha(__file__),control_source_sha256=c['source_sha256'],baseline_sha256=s['baseline_sha256'],dependency_sha256={str(path.relative_to(ROOT)):sha(path) for path in dependencies}),indent=2)+'\n')
    print(str(REPORT),flush=True)

if __name__=='__main__':main()
