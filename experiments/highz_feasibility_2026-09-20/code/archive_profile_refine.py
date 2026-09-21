#!/usr/bin/env python3
"""Recorded adaptive optimization follow-up; no residual-based pixel removal."""
import json
import numpy as np
from archive_profile_pilot import Model,fit,OUT,CONFIG
# J004131 has a weak broad feature in FeII2382 near -10 km/s, where all
# initial generic seeds had omitted gas components. Add an ordinary shared
# gas seed there before declaring the conventional model inadequate.
target='J004131-493611'
for n in [3,4,6]:
 m=Model(target,n);p,lo,hi=m.initial(0);p[n:2*n]=np.r_[-10.,np.linspace(23,49,n-1)]
 p[:n]=np.r_[12.5,np.full(n-1,14.4-np.log10(n-1))];p[2*n:3*n]=np.log(np.r_[8.,np.full(n-1,4.)])
 seed=dict(labels=m.labels,parameters=p.tolist())
 h0=fit(target,n,False,3,seed,'_blue_start',600)
 h1=fit(target,n,True,3,h0,'_blue_start',600)
 fit(target,n,False,3,h1,'_blue_cross',600)
# Regenerate comparisons from every recorded start, retaining full provenance.
for target in CONFIG:
 rows=[]
 for n in [1,2,3,4,6]:
  group=[json.loads(p.read_text()) for p in OUT.glob(f'{target}_n{n}_*.json') if '_os33' not in p.name]
  h0=min([x for x in group if not x['free_shifts']],key=lambda x:x['chi2']);h1=min([x for x in group if x['free_shifts']],key=lambda x:x['chi2'])
  rows.append(dict(ncomp=n,null=h0['name'],alternative=h1['name'],chi2_null=h0['chi2'],chi2_alternative=h1['chi2'],delta_chi2=h0['chi2']-h1['chi2'],ndata=h0['ndata'],nominal_ndf_null=h0['nominal_ndf'],aicc_null=h0['aicc'],shifts_m_s=h1['shifts_m_s'],both_success=h0['optimizer_success'] and h1['optimizer_success'],null_active_bounds=h0['active_bounds'],alternative_active_bounds=h1['active_bounds']))
 (OUT/f'{target}_comparison.json').write_text(json.dumps(dict(configuration=CONFIG[target],comparisons=rows,adaptive_followup='For J004131 counts3/4/6, an added ordinary shared component start at -10 km/s tests a weak Fe2382 feature. No line-specific opacity or pixel clipping.',interpretation='Exploratory conventional-profile feasibility; no alpha or cosmic-time-law measurement.'),indent=2)+'\n')
 best=min(rows,key=lambda x:x['aicc_null']);n=best['ncomp'];h0=json.loads((OUT/f'{best["null"]}.json').read_text());h1=json.loads((OUT/f'{best["alternative"]}.json').read_text())
 r0=fit(target,n,False,h0['seed'],h0,'_os33',600,33);r1=fit(target,n,True,h1['seed'],h1,'_os33',600,33)
 r0c=fit(target,n,False,h1['seed'],r1,'_os33_cross',600,33)
 r0=min([r0,r0c],key=lambda x:x['chi2']);r1c=fit(target,n,True,r0['seed'],r0,'_os33_cross',600,33);r1=min([r1,r1c],key=lambda x:x['chi2'])
 (OUT/f'{target}_selected_refined.json').write_text(json.dumps(dict(target=target,configuration=CONFIG[target],selection='Minimum exploratory diagonal-error null AICc among n=1,2,3,4,6; higher quadrature refinement; not a blinded model-selection protocol.',ncomp=n,null=r0['name'],alternative=r1['name'],chi2_null=r0['chi2'],chi2_alternative=r1['chi2'],delta_chi2=r0['chi2']-r1['chi2'],ndata=r0['ndata'],nominal_ndf_null=r0['nominal_ndf'],extra_shifts=r1['npar']-r0['npar'],shifts_m_s=r1['shifts_m_s'],null_active_bounds=r0['active_bounds'],alternative_active_bounds=r1['active_bounds'],both_success=r0['optimizer_success'] and r1['optimizer_success']),indent=2)+'\n')
