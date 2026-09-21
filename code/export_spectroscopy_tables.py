"""Export manuscript control comparisons from saved fits; never rerun a fit."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / 'experiments/highz_feasibility_2026-09-20'
PAPER = ROOT / 'paper_spectroscopy'
sources = {}


def read(relative):
    path = EXP / relative
    sources[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    return json.loads(path.read_text())


def controls():
    baseline = read('results/espresso_null/comparison.json')['comparisons']
    rows = []
    for case, label in [
        ('primary', 'Statistical errors, diagonal, 45 components'),
        ('expected_fluctuation', 'Expected fluctuations, diagonal, 45 components'),
    ]:
        s = next(s for s in baseline if s['case'] == case)
        rows.append((label, s['n_pixels'], s['null_parameters'], s['null_chi2'], s['delta_chi2']))
    for folder, case, label in [
        ('empirical_gls', '', 'Full pixels, GLS, 45 components'),
        ('joint_controls', 'core01_gls_', r'Joint $F<0.1$ mask + GLS'),
        ('joint_controls', 'core03_gls_', r'Joint padded $F<0.3$ mask + GLS'),
        ('gas_structure', '', 'Full pixels, GLS, 40 components'),
    ]:
        s = read(f'results/{folder}/{case}comparison.json')
        n = read(f"results/{folder}/{s['null']}.json")
        rows.append((label, n['ndata'], n['npar'], n['chi2'], s['delta_chi2']))
    lines = [r'\begin{tabular}{lrrrr}', r'\toprule',
             r'Selection and gas architecture & Pixels & $k_0$ & $Q_0$ & $\Delta Q$ \\',
             r'\midrule']
    for label, pixels, k, q0, dq in rows:
        lines.append(f'{label} & {pixels} & {k} & {q0:.2f} & {dq:.2f} ' + r'\\')
    lines += [r'\bottomrule', r'\end{tabular}', '']
    (PAPER / 'control_table.tex').write_text('\n'.join(lines))
    return rows


def archive():
    rows = sorted(read('results/archive_complete/six_target_comparison.json')['rows'],
                  key=lambda row: row['z_abs'])
    lines = [r'\begingroup', r'\setlength{\tabcolsep}{4pt}',
             r'\begin{tabular}{lrrrrrl}', r'\toprule',
             r'Target ($z$) & Pixels & $n_g$ & $Q_0$ & $\Delta Q$ & $k_\delta$ & Status \\',
             r'\midrule']
    for r in rows:
        target = r['target'].replace('-', '$-$')
        status = 'Conditional' if r['adequacy_flag'] else 'Inadequate'
        lines.append(f"{target} ({r['z_abs']:.3f}) & {r['pixels']} & {r['ncomp']} & "
                     f"{r['Q0']:.2f} & {r['delta_Q']:.2f} & {r['extra_shifts']} & {status} " + r'\\')
    lines += [r'\bottomrule', r'\end{tabular}', r'\endgroup', '']
    (PAPER / 'archive_table.tex').write_text('\n'.join(lines))
    return rows


if __name__ == '__main__':
    rows = controls()
    archive_rows = archive()
    (PAPER / 'table_provenance.json').write_text(json.dumps({
        'purpose': 'Saved comparison export; no new optimization or significance calculation',
        'control_rows': rows,
        'archive_rows': archive_rows,
        'experiment_source_sha256': sources,
        'export_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }, indent=2) + '\n')
    print(f'Exported {len(rows)} saved control comparisons and {len(archive_rows)} archive pairs.')
