#!/usr/bin/env python3
"""Export the actual-data table and figures into the parent manuscript."""
from pathlib import Path
import csv
import shutil

BASE=Path(__file__).resolve().parents[1]
ROOT=BASE.parents[1]


def run():
    rows=list(csv.DictReader((BASE/'results/observed_profile_metrics.csv').open()))
    lines=[r'\begin{tabular}{llrrrr}',r'\toprule',
        r'Sightline & Fe II & $W_r$ (m\AA) & $\bar v$ (km/s) & $\Delta_C$ (km/s) & $\Delta_W$ (km/s) \\',r'\midrule']
    for r in rows:
        if r['status']!='descriptive_profile_only':continue
        target='HE 0515' if r['target']=='J051707-441055' else 'Q0347'
        ion=r['transition'].split()[-1]
        w=float(r['rest_EW_mA']);ws=float(r['rest_EW_error_diagonal_mA'])
        v=float(r['centroid_km_s']);vs=float(r['centroid_error_diagonal_km_s'])
        dc=float(r['centroid_continuum_1pct_max_change_km_s']);dw=float(r['centroid_range_10kms_max_change_km_s'])
        lines.append(f'{target} & {ion} & ${w:.2f}\\pm{ws:.2f}$ & ${v:.3f}\\pm{vs:.3f}$ & {dc:.3f} & {dw:.3f} '+r'\\')
    lines += [r'\bottomrule',r'\end{tabular}']
    (ROOT/'paper/highz_observed_table.tex').write_text('\n'.join(lines)+'\n')
    for source,dest in [('observed_profile_diagnostics','highz_observed_profiles'),
                        ('observed_registration','highz_observed_registration')]:
        for ext in ('pdf','png'):
            shutil.copyfile(BASE/f'figures/{source}.{ext}',ROOT/f'figures/{dest}.{ext}')
    print(ROOT/'paper/highz_observed_table.tex')


if __name__=='__main__':run()
