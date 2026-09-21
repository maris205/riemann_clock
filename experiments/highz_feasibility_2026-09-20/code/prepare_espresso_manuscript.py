#!/usr/bin/env python3
"""Export final selected comparison without changing third-party inputs."""
from pathlib import Path
import json,shutil
from summarize_espresso_null import main,select,ROOT,OUT
PAPER=ROOT.parents[1]/'paper';FIG=ROOT.parents[1]/'figures'
NAMES={'primary':'Primary errors','expected_fluctuation':'Expected fluctuation','ar1_rho03':r'Assumed AR(1), $\rho=0.3$','lsf_minus3pct':r'LSF width $-3\%$','lsf_plus3pct':r'LSF width $+3\%$','omit_blended':'Omit 2344, 2586','opacity_zero':'Opacity and zero freedom','omit_saturated':'Omit 2382, 2600','wider_bounds':'Wider nuisance bounds'}
def export():
 main();x=json.loads((OUT/'comparison.json').read_text());a=select('null');b=select('alternative')
 macros={'EspressoNullChi':f"{a['chi2']:.2f}",'EspressoNullNu':str(a['nominal_ndf']),'EspressoDeltaChi':f"{a['chi2']-b['chi2']:.2f}",'EspressoAdded':str(b['npar']-a['npar'])}
 (PAPER/'espresso_null_numbers.tex').write_text('% Generated from selected fitted results; do not edit numbers manually.\n'+''.join('\\newcommand{\\'+k+'}{'+v+'}\n' for k,v in macros.items()))
 lines=[r'\begin{tabular}{lrrrr}',r'\toprule',r'Case & Pixels & $\chi^2_0/\nu_{\rm nom}$ & $\Delta\chi^2$ & Extra shifts\\',r'\midrule']
 for r in x['comparisons']:
  if not r['null_success'] or not r['alternative_success']:raise RuntimeError('An unconverged selected fit must not enter the final comparison table')
  lines.append(f"{NAMES[r['case']]} & {r['n_pixels']} & {r['null_chi2_per_nominal_ndf']:.3f} & {r['delta_chi2']:.2f} & {r['added_parameters']}"+r'\\')
 lines += [r'\bottomrule',r'\end{tabular}']
 (PAPER/'espresso_null_table.tex').write_text('\n'.join(lines)+'\n')
 for n in ['espresso_conventional_null','espresso_null_robustness']:
  for ext in ['pdf','png']:shutil.copy2(ROOT/'figures'/f'{n}.{ext}',FIG/f'{n}.{ext}')
 (OUT/'selected_primary.json').write_text(json.dumps({'null':a['name'],'alternative':b['name'],'chi2_reduction':a['chi2']-b['chi2'],'shifts_m_s':{k:1000*v for k,v in b['shifts_km_s'].items()}},indent=2)+'\n')
if __name__=='__main__':export()
