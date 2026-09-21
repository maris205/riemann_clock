#!/usr/bin/env python3
"""Export completed nonlinear controls; never run or overwrite optimizations."""
from pathlib import Path
import json
import hashlib
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
PAPER = PROJECT / 'paper'
FIG = PROJECT / 'figures'
OUT = ROOT / 'results/followup_summary'


def read(path):
    return json.loads(path.read_text())


def pair(folder, summary_name, case, label):
    directory = ROOT / 'results' / folder
    s = read(directory / summary_name)
    a = read(directory / (s['null'] + '.json'))
    b = read(directory / (s['alternative'] + '.json'))
    assert a['ndata'] == b['ndata']
    assert a['optimizer_success'] and b['optimizer_success']
    assert abs(a['chi2'] - b['chi2'] - s['delta_chi2']) < 1e-8
    return dict(case=case, label=label, null_file=str((directory / (s['null'] + '.json')).relative_to(ROOT)),
                alternative_file=str((directory / (s['alternative'] + '.json')).relative_to(ROOT)),
                ndata=a['ndata'], added_parameters=b['npar'] - a['npar'],
                null_chi2=a['chi2'], alternative_chi2=b['chi2'], delta_chi2=s['delta_chi2'],
                shifts_m_s=b['shifts_m_s'], null_optimality=a['optimality'],
                alternative_optimality=b['optimality'], null_stop=a['message'], alternative_stop=b['message'])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    old = read(ROOT / 'results/espresso_null/comparison.json')['comparisons']
    original = next(r for r in old if r['case'] == 'primary')
    primary = dict(case='primary', label='Original diagonal errors', ndata=original['n_pixels'],
                   added_parameters=original['added_parameters'], null_chi2=original['null_chi2'],
                   alternative_chi2=original['alternative_chi2'], delta_chi2=original['delta_chi2'])
    gls = pair('empirical_gls', 'comparison.json', 'gls', 'Continuum-informed GLS')
    core1 = pair('saturation_controls', 'core01_comparison.json', 'core01', 'Mask strong-line F < 0.1')
    core3 = pair('saturation_controls', 'core03_padded_comparison.json', 'core03', 'Mask F < 0.3 + 3 pixels')
    rows = [primary, gls, core1, core3]
    control = read(ROOT / 'results/noise_covariance/controls.json')['primary']
    summary = dict(espresso_comparisons=rows, continuum=control,
                   interpretation='Separate exploratory controls, not their joint application. Differences are conditional nonlinear objective reductions, not calibrated physical-discovery significances. Optimizer stopping criteria do not establish stationarity or global minima.')
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    macros = {'EspressoGLSDelta': gls['delta_chi2'], 'EspressoCoreOneDelta': core1['delta_chi2'],
              'EspressoCoreThreeDelta': core3['delta_chi2']}
    (PAPER / 'espresso_followup_numbers.tex').write_text(
        '% Generated from completed paired fits; see prepare_followup_manuscript.py.\n' +
        ''.join('\\newcommand{\\' + k + '}{' + f'{v:.2f}' + '}\n' for k, v in macros.items()))
    labels = ['Original diagonal errors', 'Continuum-informed GLS',
              r'Mask strong-line $F<0.1$', r'Mask $F<0.3$, pad 3 pixels']
    lines = [r'\begin{tabular}{lrrrr}', r'\hline',
             r'ESPRESSO comparison & Pixels & $\chi^2_0$ & $\chi^2_1$ & $\Delta\chi^2$ \\', r'\hline']
    for r, label in zip(rows, labels):
        lines.append(f'{label} & {r["ndata"]} & {r["null_chi2"]:.2f} & {r["alternative_chi2"]:.2f} & {r["delta_chi2"]:.2f} ' + r'\\')
    lines += [r'\hline', r'\end{tabular}']
    (PAPER / 'espresso_followup_table.tex').write_text('\n'.join(lines) + '\n')

    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.3), gridspec_kw={'width_ratios': [1., 1.4]})
    acf = np.array(control['acf']); lag = np.arange(1, len(acf))
    ci = np.array(control['interval_cluster_bootstrap']['acf_95_percentile_intervals'])
    ax = axes[0]
    ax.errorbar(lag, acf[1:], yerr=[acf[1:] - ci[1:, 0], ci[1:, 1] - acf[1:]],
                fmt='o', color='#17668b', ms=4, capsize=2, label='Continuum controls (95% interval)')
    ax.plot(lag, acf[1]**lag, '--', color='#ba752a', label='AR(1), matched lag 1')
    ax.plot(lag, acf[1:] * (1 - lag / len(acf)), ':s', color='#4f8452', ms=3, label='Fixed tapered kernel used in GLS')
    ax.axhline(0, color='.6', lw=.7)
    ax.set(xlabel='Pixel lag (0.4 km/s per pixel)', ylabel='Autocorrelation',
           title='25 disjoint continuum intervals / 13,750 pixels', xlim=(.5, 10.5))
    ax.legend(fontsize=7.6, loc='upper right')
    ax = axes[1]
    values = [r['delta_chi2'] for r in rows]
    ax.barh(np.arange(4), values, color=['#8c989d', '#17668b', '#559875', '#a78044'], height=.62)
    for i, (r, value) in enumerate(zip(rows, values)):
        ax.text(value + .65, i, f'{value:.2f}', va='center', fontsize=9)
    ax.set(yticks=np.arange(4), yticklabels=['Original errors\n2,931 pixels', 'Empirical GLS\n2,931 pixels',
                                          'Mask F < 0.1\n' + f'{core1["ndata"]:,} pixels',
                                          'Mask F < 0.3 + padding\n' + f'{core3["ndata"]:,} pixels'],
           xlabel=r'Nonlinear $\chi^2_0-\chi^2_1$ (five added shifts)',
           title='Separate paired ESPRESSO re-fits', xlim=(0, max(values) * 1.16))
    ax.invert_yaxis()
    fig.text(.5, .015, 'Conditional comparisons: masks and noise models differ; this is not a discovery significance or a joint-controls test.', ha='center', fontsize=8.5)
    fig.tight_layout(rect=[0, .06, 1, 1])
    for extension in ['pdf', 'png']:
        fig.savefig(OUT / f'espresso_noise_controls.{extension}', dpi=180, bbox_inches='tight')
        shutil.copyfile(OUT / f'espresso_noise_controls.{extension}', FIG / f'espresso_noise_controls.{extension}')
    plt.close(fig)
    dependencies = [ROOT / 'results/noise_covariance/controls.json']
    for r in rows[1:]:
        dependencies.extend(ROOT / r[k] for k in ['null_file', 'alternative_file'])
    manifest = {str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies}
    (OUT / 'source_hashes.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
