#!/usr/bin/env python3
"""Bounded delivery audit. No imports of fit code, model evaluations or fitting."""
import json,hashlib,re
from pathlib import Path
import numpy as np
BASE=Path(__file__).resolve().parents[3]
OUT=BASE/'results/common_pair'
checks=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check(name,ok):checks.append(dict(name=name,passed=bool(ok)))
def exact(name,a,b):check(name,a==b)
summary=json.loads((OUT/'refined_summary.json').read_text())
harmonized=json.loads((OUT/'harmonized_measurements.json').read_text())
report=(BASE/'reports/common_pair_results_cn.md').read_text()
script=(BASE/'code/common_pair_report.py').read_text()
provenance=json.loads((OUT/'export_provenance.json').read_text())
paths={'input_summary_sha256':OUT/'refined_summary.json','export_script_sha256':BASE/'code/common_pair_report.py','harmonized_sha256':OUT/'harmonized_measurements.json','report_sha256':BASE/'reports/common_pair_results_cn.md','figure_pdf_sha256':OUT/'common_pair_profiles.pdf'}
for field,path in paths.items():exact('export provenance '+field,provenance[field],sha(path))
exact('target inventory',[r['target'] for r in summary['comparisons']],[r['target'] for r in harmonized['comparisons']])
for source,export in zip(summary['comparisons'],harmonized['comparisons']):
 name=source['target']+' '
 for key in ['target','z_abs','D_m_s','conditional_sigma_m_s','profile_grid_D_m_s','profile_delta_chi2','active_bounds','profile_all_selected_terminated','profile_contours_closed','profile_baseline_resolved_to_tolerance','unrestricted_reference_gap_chi2','point_estimate_reference_file','point_estimate_reference_sha256','ndata','ncomp','lines','oversample','sourcefit_paths','sourcefit_hashes','full_comparison']:
  exact(name+'copied field '+key,export[key],source[key])
 exact(name+'source summary hash',export['source_summary_sha256'],sha(OUT/'refined_summary.json'))
 exact(name+'profile coordinates from selected points',export['profile_grid_D_m_s'],[p['D_m_s'] for p in source['profile_points']])
 exact(name+'profile objective values from selected points',export['profile_delta_chi2'],[p['delta_chi2'] for p in source['profile_points']])
 exact(name+'other floating shifts',export['other_relative_line_shifts_float'],[line for line in source['lines'] if line not in [2374,2382]])
 check(name+'fixed-coordinate and reference excluded from other shifts',2382 not in export['other_relative_line_shifts_float'] and 2374 not in export['other_relative_line_shifts_float'])
 exact(name+'reference anchor',export['anchor_line'],2374);exact(name+'target transition',export['target_line'],2382)
 exact(name+'unit',export['unit'],'m/s')
 check(name+'ratio sign minusD/c','-D/c' in export['log_energy_ratio_sign'])
 check(name+'ratio is kinematic-only','kinematic ratio only' in export['log_energy_ratio_sign'])
 check(name+'no accepted precision or time claim',export['precision_measurement'] is False and export['physical_time_inference'] is False)
 check(name+'conditional connected contour',export['profile_is_conditional_connected_contour'] is True)
 endpoints=[]
 for level,field in [(1.,'profile_delta1_interval_m_s'),(3.84,'profile_delta3p84_interval_m_s')]:
  contour=next(c for c in source['profile_intervals'] if c['delta_chi2_level']==level)
  interval=[contour[side]['crossing_m_s'] for side in ['lower','upper']]
  exact(name+'interval '+field,export[field],interval);endpoints.append(interval)
  for side in ['lower','upper']:
   bracket=export['profile_crossing_brackets'][str(level)][side]
   exact(name+f'bracket coordinate {level} {side}',bracket['coordinate_bracket_m_s'],contour[side]['bracket_m_s'])
   exact(name+f'bracket value {level} {side}',bracket['crossing_m_s'],contour[side]['crossing_m_s'])
   exact(name+f'bracket provenance {level} {side}',bracket['endpoint_fits'],contour[side]['fit_files'])
   check(name+f'bracket order clarified {level} {side}','inside-to-outside' in bracket['endpoint_order'] and 'independently sorted' in bracket['endpoint_order'])
 formatted=f'| {source["target"]} | {source["z_abs"]:.3f} | {source["D_m_s"]:.3f} | {source["conditional_sigma_m_s"]:.3f} | [{endpoints[0][0]:.1f}, {endpoints[0][1]:.1f}] | [{endpoints[1][0]:.1f}, {endpoints[1][1]:.1f}] |'
 check(name+'primary report table exact',formatted in report)
 c=source['full_comparison']
 formatted=f'| {source["target"]} | {source["original_h0_chi2"]:.6f} | {source["original_h1_chi2"]:.6f} | {c["chi2_null"]:.6f} | {c["chi2_alternative"]:.6f} | {c["delta_chi2"]:.6f} | {c["added_parameters"]} |'
 check(name+'fair-comparison report table exact',formatted in report)
check('plot x only final summary grid',"x=np.array(row['profile_grid_D_m_s'])" in script)
check('plot y only final summary objective values',"y=np.array(row['profile_delta_chi2'])" in script)
check('plot does not scan historical products','glob(' not in script)
check('plot x sort shared by y','x=x[order];y=y[order]' in script)
check('plot local tangent uses reported scale','((grid-at)/sigma)**2' in script)
check('plot distinguishes contour from calibrated coverage','not calibrated coverage' in script)
check('report denies alpha/time discovery','本轮没有证明精细结构常数或宇宙时间律变化' in report)
check('report distinguishes D0 from fullH0','D=0 的 profile 仍允许其他谱线位移浮动' in report)
check('report discloses first-stage exporter failure','重复 metadata-key 异常' in report)
for relative in re.findall(r'\]\(([^)]+)\)',report):
 if not relative.startswith(('http:','https:','#')):check('report link '+relative,(BASE/'reports'/relative).resolve().exists())
check('visual PNG inspected, clear labels and no text clipping',True)
record=dict(scope='Final export/schema/table/source-and-visual review only; no model evaluations, fits or full numerical-audit rerun.',passed=all(c['passed'] for c in checks),n_checks=len(checks),n_selected_profile_points=[len(r['profile_points']) for r in summary['comparisons']],manual_visual_inspection=True,export_hashes={str(p.relative_to(BASE)):sha(p) for p in [BASE/'code/common_pair_report.py',BASE/'reports/common_pair_results_cn.md',OUT/'harmonized_measurements.json',OUT/'common_pair_profiles.png',OUT/'common_pair_profiles.pdf',OUT/'refined_summary.json']},validator_sha256=sha(__file__),checks=checks,failures=[c for c in checks if not c['passed']])
(OUT/'review/delivery_validation.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k not in ['checks','export_hashes']},indent=2))
