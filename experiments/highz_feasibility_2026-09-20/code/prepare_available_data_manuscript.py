#!/usr/bin/env python3
"""Export the completed combined-noise/mask and gas-architecture comparisons."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
PAPER = PROJECT/'paper'
FIG = PROJECT/'figures'
OUT = ROOT/'results/available_data_campaign'


def read(path):
    return json.loads(path.read_text())


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sources = []
    def source(path):
        sources.append(path)
        return read(path)
    full = source(ROOT/'results/empirical_gls/comparison.json')
    old = source(ROOT/'results/espresso_null/comparison.json')['comparisons']
    baseline = next(r for r in old if r['case']=='primary')
    mild = source(ROOT/'results/joint_controls/core01_gls_comparison.json')
    wide = source(ROOT/'results/joint_controls/core03_gls_comparison.json')
    gas = source(ROOT/'results/gas_structure/comparison.json')
    prior = [source(ROOT/'results/saturation_controls'/f'{case}_comparison.json')
             for case in ['core01','core03_padded']]
    info = source(ROOT/'results/mask_information/information.json')
    infos = [next(x for x in c['cutoffs'] if x['relative_svd_cutoff']==1e-10) for c in info['cases']]
    retention = [x['coherent_information_retention_vs_full'] for x in infos]
    rows = []
    configs = [('Full pixels, 45 components', full, 'empirical_gls', 45),
               (r'Joint $F<0.1$ mask + GLS', mild, 'joint_controls', 45),
               (r'Joint padded $F<0.3$ mask + GLS', wide, 'joint_controls', 45),
               ('Full pixels, 40 components', gas, 'gas_structure', 40)]
    for title, s, directory, ncomp in configs:
        a=source(ROOT/'results'/directory/(s['null']+'.json'))
        b=source(ROOT/'results'/directory/(s['alternative']+'.json'))
        assert a['ndata']==b['ndata'] and b['npar']-a['npar']==5
        assert abs(a['chi2']-b['chi2']-s['delta_chi2'])<1e-8
        rows.append(dict(label=title,ncomp=ncomp,ndata=a['ndata'],null_parameters=a['npar'],
                         alternative_parameters=b['npar'],null_chi2=a['chi2'],alternative_chi2=b['chi2'],
                         delta_chi2=s['delta_chi2'],null_success=a['optimizer_success'],
                         alternative_success=b['optimizer_success']))
    (OUT/'model_comparisons.json').write_text(json.dumps(rows,indent=2)+'\n')
    table=[r'\begin{tabular}{lrrrr}',r'\hline',
           r'GLS comparison & Pixels & $k_0$ & $\chi^2_0$ & $\Delta\chi^2$ \\',r'\hline']
    for row in rows:
        table.append(f'{row["label"]} & {row["ndata"]} & {row["null_parameters"]} & {row["null_chi2"]:.2f} & {row["delta_chi2"]:.2f} '+r'\\')
    table += [r'\hline',r'\end{tabular}']
    (PAPER/'available_data_table.tex').write_text('\n'.join(table)+'\n')
    macros=dict(JointCoreOneDelta=mild['delta_chi2'],JointCoreThreeDelta=wide['delta_chi2'],GasFortyDelta=gas['delta_chi2'])
    (PAPER/'available_data_numbers.tex').write_text('% Generated from completed saved comparisons.\n'+''.join(
        '\\newcommand{\\'+key+'}{'+f'{val:.2f}'+'}\n' for key,val in macros.items()))
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),gridspec_kw={'width_ratios':[1.25,1]})
    x=np.arange(3)
    diagonal=[baseline['delta_chi2']]+[a['delta_chi2'] for a in prior]
    gls=[full['delta_chi2'],mild['delta_chi2'],wide['delta_chi2']]
    ax=axes[0]
    for pos,values,color,label in [(-.18,diagonal,'#939caa','Diagonal errors'),(.18,gls,'#216b8d','Empirical covariance')]:
        bars=ax.bar(x+pos,values,width=.34,color=color,label=label)
        ax.bar_label(bars,labels=[f'{v:.2f}' for v in values],fontsize=9,padding=3)
    labels=['Full pixels\n2,931','Dark-core mask\n2,692','Wider padded mask\n2,473']
    ax.set(xticks=x,xticklabels=labels,ylabel=r'Nonlinear $\chi^2_0-\chi^2_1$',ylim=(0,61),title='Actual paired fits: five extra shifts')
    ax.legend(fontsize=9)
    ax=axes[1]
    bars=ax.bar(x,np.array(retention)*100,color=['#216b8d','#529878','#ad884a'],width=.6)
    ax.bar_label(bars,labels=[f'{100*v:.1f}%' for v in retention],padding=4)
    ax.set(xticks=x,xticklabels=['Full GLS','Dark-core\nmask + GLS','Wider mask\n+ GLS'],ylim=(0,116),
           ylabel='Local information retained (%)',title='Illustrative coherent strong-line direction')
    ax.text(.98,.87,'2382 and 2600 each +100 m/s\nCommon fixed gas parameters\nLinear nuisance projection',
            ha='right',va='top',transform=ax.transAxes,fontsize=8.5)
    fig.text(.5,.015,'Left: observed conditional objective reductions. Right: a local sensitivity calculation, not an observed detection statistic.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=[0,.055,1,1])
    for extension in ['pdf','png']:
        fig.savefig(FIG/f'joint_controls_information.{extension}',dpi=180,bbox_inches='tight')
    plt.close(fig)
    (OUT/'model_source_hashes.json').write_text(json.dumps({str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                                         for p in dict.fromkeys(sources)},indent=2)+'\n')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    main()
