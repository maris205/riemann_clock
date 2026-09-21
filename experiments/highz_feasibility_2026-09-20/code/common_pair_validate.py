#!/usr/bin/env python3
"""Independent common Fe II pair numerical audit; performs no optimization."""
from pathlib import Path
import argparse
import hashlib
import json
import tempfile
from types import SimpleNamespace
import numpy as np
from scipy.linalg import qr
import archive_profile_pilot as pilot
import archive_profile_expected_noise as expected

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/common_pair/review'
CHECKS = []


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check(name, passed, detail=None):
    CHECKS.append(dict(name=name, passed=bool(passed), detail=detail))


def close(name, actual, reference, rtol=2e-9, atol=2e-10):
    a, b = np.asarray(actual), np.asarray(reference)
    ok = a.shape == b.shape and np.allclose(a, b, rtol=rtol, atol=atol)
    difference = float(np.max(np.abs(a.astype(float)-b.astype(float)))) if a.shape == b.shape and a.size else None
    check(name, ok, difference)


def curvature(J, index):
    """Independent QR projection and full normalized-J SVD inverse covariance."""
    nuisance = np.delete(J, index, axis=1)
    ns = np.linalg.norm(nuisance, axis=0)
    Q, R, permutation = qr(nuisance/ns, mode='economic', pivoting=True)
    # These two problems have full column rank at the validated cutoff.
    sv_n = np.linalg.svd(nuisance/ns, compute_uv=False)
    rank = int(np.sum(sv_n > sv_n[0]*1e-10))
    q = J[:, index] - Q[:, :rank] @ (Q[:, :rank].T @ J[:, index])
    qr_sigma = 1000 / np.linalg.norm(q)
    scale = np.linalg.norm(J, axis=0)
    U, s, Vh = np.linalg.svd(J/scale, full_matrices=False)
    full_rank = int(np.sum(s > s[0]*1e-10))
    covariance_diagonal = np.sum((Vh[:full_rank, index]/s[:full_rank])**2)/scale[index]**2
    svd_sigma = 1000*np.sqrt(covariance_diagonal)
    return dict(sigma_qr_m_s=float(qr_sigma), sigma_svd_m_s=float(svd_sigma),
                nuisance_rank=rank, full_rank=full_rank,
                normalized_singular_ratio=float(s[-1]/s[0]))


def baseline_review():
    original = ROOT/'results/archive_expansion/profile_pilot/expected_noise'
    summary = json.loads((original/'summary.json').read_text())
    snapshot = json.loads((OUT/'baseline_snapshot.json').read_text())
    for filename, digest in snapshot.items():
        check('unchanged baseline '+filename, sha(ROOT/filename) == digest)
    rows = []
    for row in summary['comparisons']:
        path = original/(row['expected_alternative']+'.json')
        record = json.loads(path.read_text())
        saved = np.load(path.with_suffix('.npz'))
        p = np.array(record['parameters'])
        m = expected.ExpectedModel(record['target'],record['ncomp'],True,33)
        prefix = record['target']+' original endpoint '
        r, J, _ = m.evaluate(p)
        close(prefix+'replayed residual', r, saved['residual'])
        close(prefix+'replayed Jacobian', J, saved['jacobian'])
        close(prefix+'replayed objective',r@r,record['chi2'])
        check(prefix+'fixed 2374 reference',2374 not in m.shiftkeys and 'shift_2374' not in m.labels)
        check(prefix+'all remaining H1 shifts present',set(m.shiftkeys) == set(m.keys)-{2374})
        index = m.labels.index('shift_2382')
        info = curvature(J,index)
        check(prefix+'full nuisance rank',info['nuisance_rank']==m.npar-1)
        check(prefix+'full joint rank',info['full_rank']==m.npar)
        close(prefix+'QR versus full-SVD sigma',info['sigma_qr_m_s'],info['sigma_svd_m_s'],rtol=1e-11)
        D = 1000*p[index]
        close(prefix+'declared contrast',D,record['shifts_m_s']['2382'])
        # Independent central differences include the target and every type of nuisance.
        names = ['shift_2382','logN_0','v_0','logb_0','cont0_2374',
                 'cont1_2382','zero_2382','logfwhm_2382',f'shift_{m.shiftkeys[0]}']
        finite_differences = []
        for name in names:
            j = m.labels.index(name)
            step = 2e-5 if name.startswith(('v_','shift_')) else 2e-6
            a=p.copy(); b=p.copy(); a[j]+=step; b[j]-=step
            fd=(m.fun(a)-m.fun(b))/(2*step)
            rel=float(np.linalg.norm(fd-J[:,j])/max(np.linalg.norm(fd),1e-10))
            # logfwhm is itself a centred numerical convolution derivative.
            tolerance=3e-4 if name.startswith('logfwhm') else 3e-5
            check(prefix+'finite difference '+name,rel<tolerance,rel)
            finite_differences.append(dict(label=name,relative_norm_error=rel))
        # The formula must use supplied absolute row2 noise, not rescale by reduced chi².
        rows.append(dict(target=m.target,contrast_m_s=float(D),**info,
                         active_nuisance_bounds=record['active_bounds'],
                         finite_differences=finite_differences))
    return rows


def optimized_gradient(J, residual, p, lo, hi, indices):
    gradient=J[:,indices].T@residual
    position=np.asarray(p)[indices]
    low=np.asarray(lo)[indices]; high=np.asarray(hi)[indices]
    distance=np.ones(len(indices))
    distance[gradient<0]=(high-position)[gradient<0]
    distance[gradient>0]=(position-low)[gradient>0]
    return float(np.max(np.abs(gradient*distance)))


def actual_wrapper_probes():
    """Execute the actual optimizer closures with a recorder, not an optimizer."""
    import common_pair_fit as fitcode
    summary=json.loads((expected.OUT/'summary.json').read_text())
    original_out=fitcode.OUT
    original_optimizer=fitcode.least_squares
    probes=[]
    try:
        with tempfile.TemporaryDirectory(prefix='wrapper_probe_',dir=OUT) as temporary:
            fitcode.OUT=Path(temporary)
            for selection in summary['comparisons']:
                model=expected.ExpectedModel(selection['target'],selection['ncomp'],True,33)
                original=json.loads((expected.OUT/(selection['expected_alternative']+'.json')).read_text())
                at=model.labels.index('shift_2382')
                for fixed in [None,0.,-1.,1.]:
                    prefix=selection['target']+' actual wrapper '+str(fixed)+' '
                    expected_free=np.array([j for j in range(model.npar) if fixed is None or j!=at])
                    p,lo,hi=model.initial()
                    mapping=dict(zip(original['labels'],original['parameters']))
                    p=np.clip([mapping[k] for k in model.labels],lo+1e-9,hi-1e-9)
                    if fixed is not None:p[at]=fixed
                    invoked=[]
                    def recorder(fun,x0,jac,bounds,**options):
                        invoked.append(True)
                        close(prefix+'initial optimizer vector',x0,p[expected_free],rtol=0,atol=0)
                        close(prefix+'optimizer lower bounds',bounds[0],lo[expected_free],rtol=0,atol=0)
                        close(prefix+'optimizer upper bounds',bounds[1],hi[expected_free],rtol=0,atol=0)
                        points=[np.asarray(x0).copy()]
                        test=np.asarray(x0).copy()
                        direction=np.where(test<(lo[expected_free]+hi[expected_free])/2,1.,-1.)
                        test+=direction*1e-6
                        points.append(test)
                        for k,x in enumerate(points):
                            independent=p.copy();independent[expected_free]=x
                            r,J,_=model.evaluate(independent)
                            close(prefix+'closure residual '+str(k),fun(x),r)
                            close(prefix+'closure Jacobian '+str(k),jac(x),J[:,expected_free])
                        r,J,_=model.evaluate(p)
                        optimality=optimized_gradient(J,r,p,lo,hi,expected_free)
                        return SimpleNamespace(x=x0,success=False,status=0,message='Independent closure probe; no optimization',nfev=0,optimality=optimality)
                    fitcode.least_squares=recorder
                    name='probe_'+('free' if fixed is None else str(fixed).replace('-','minus').replace('.','_'))
                    record=fitcode.fit(name,model,expected.OUT/(selection['expected_alternative']+'.json'),fixed,
                                       fitcode.metadata(model,selection),max_nfev=1)
                    check(prefix+'recorder invoked exactly once',len(invoked)==1)
                    close(prefix+'stored full parameters',record['parameters'],p,rtol=0,atol=0)
                    close(prefix+'stored optimizer indices',record['optimized_parameter_indices'],expected_free,rtol=0,atol=0)
                    probes.append(dict(target=model.target,fixed_shift_km_s=fixed,npar_optimized=len(expected_free)))
    finally:
        fitcode.OUT=original_out
        fitcode.least_squares=original_optimizer
    return probes


def fit_products_review():
    paths=sorted((ROOT/'results/common_pair').glob('J*/*.json'))
    records={}
    models={}
    for path in paths:
        obj=json.loads(path.read_text())
        if 'optimized_parameter_indices' not in obj:continue
        target=obj['target'];prefix=target+' '+obj['name']+' '
        free='shift_2382' in obj['labels'];model_key=(target,free)
        if model_key not in models:models[model_key]=expected.ExpectedModel(target,obj['ncomp'],free,33)
        model=models[model_key]
        p=np.asarray(obj['parameters']);r,J,profiles=model.evaluate(p)
        saved=np.load(path.with_suffix('.npz'))
        rel=str(path.relative_to(ROOT));records[rel]=obj
        at=model.labels.index('shift_2382') if free else None;fixed=obj['fixed_shift_km_s']
        indices=np.array([i for i in range(model.npar) if fixed is None or i!=at])
        _,lo,hi=model.initial()
        check(prefix+'labels and dimension',obj['labels']==model.labels and obj['ndata']==model.ndata and obj['npar_total']==model.npar)
        close(prefix+'optimized index mapping',obj['optimized_parameter_indices'],indices,rtol=0,atol=0)
        close(prefix+'saved optimized index mapping',saved['optimized_parameter_indices'],indices,rtol=0,atol=0)
        check(prefix+'optimized dimension',obj['npar_optimized']==len(indices))
        if fixed is not None:close(prefix+'fixed coordinate',p[at],fixed,rtol=0,atol=0)
        close(prefix+'common contrast',obj['D_m_s'],1000*p[at] if free else 0,rtol=0,atol=0)
        close(prefix+'lower bounds unchanged',obj['bounds_lower'],lo,rtol=0,atol=0)
        close(prefix+'upper bounds unchanged',obj['bounds_upper'],hi,rtol=0,atol=0)
        check(prefix+'final feasible',np.all(p>=lo) and np.all(p<=hi))
        start_path=ROOT/obj['start_file'];start=json.loads(start_path.read_text())
        lookup=dict(zip(start['labels'],start['parameters']));default,_,_=model.initial()
        initial=np.clip([lookup.get(k,v) for k,v in zip(model.labels,default)],lo+1e-9,hi-1e-9)
        if fixed is not None:initial[at]=fixed
        close(prefix+'initial mapping',obj['initial_parameters'],initial,rtol=0,atol=0)
        check(prefix+'start hash',sha(start_path)==obj['start_sha256'])
        check(prefix+'all source hashes',all(sha(ROOT/key)==value for key,value in obj['source_sha256'].items()))
        check(prefix+'sampling 33 and fixed reference',obj['oversample']==33 and obj['anchor_line']==2374 and obj['target_line']==2382)
        close(prefix+'saved parameters',saved['parameters'],p,rtol=0,atol=0)
        close(prefix+'saved residual',saved['residual'],r)
        close(prefix+'saved Jacobian',saved['jacobian'],J)
        close(prefix+'objective',obj['chi2'],r@r)
        close(prefix+'bound-scaled gradient',obj['optimality'],optimized_gradient(J,r,p,lo,hi,indices),rtol=2e-6,atol=2e-8)
        check(prefix+'stopping flag consistent',obj['optimizer_success']==(obj['optimizer_status']>0))
        check(prefix+'evaluation budget respected',obj['nfev']<=obj['max_nfev'])
        active=[model.labels[i] for i in range(model.npar) if min(p[i]-lo[i],hi[i]-p[i])<1e-5]
        check(prefix+'active bounds',obj['active_bounds']==active)
        for line,profile in zip(model.lines,profiles):
            key=str(line['key'])
            for column in ['v','wave','flux','error','good','source_indices']:
                close(prefix+key+' '+column,saved[key+'_'+column],line[column],rtol=0,atol=0)
            close(prefix+key+' model',saved[key+'_model'],profile)
    return records


def summary_review(records,partial):
    rows=[]
    for target in pilot.CONFIG:
        folder=ROOT/'results/common_pair'/target
        refined=(folder/'refined_summary.json').exists()
        summary_path=folder/('refined_summary.json' if refined else 'summary.json')
        if not summary_path.exists():
            if not partial:check(target+' final summary exists',False)
            continue
        summary=json.loads(summary_path.read_text());prefix=target+' final summary '
        model=expected.ExpectedModel(target,summary['ncomp'],True,33)
        reference_path=ROOT/summary['point_estimate_reference_file']
        reference=json.loads(reference_path.read_text())
        full=[json.loads((ROOT/file).read_text()) for file in summary['full_fit_files']]
        fixed=[records[file] for file in summary['all_profile_fit_files']]
        all_endpoints=full+fixed
        close(prefix+'lowest of every endpoint',summary['reference_chi2'],min(x['chi2'] for x in all_endpoints),rtol=0,atol=1e-11)
        close(prefix+'selected reference objective',summary['reference_chi2'],reference['chi2'])
        check(prefix+'reference hash',sha(reference_path)==summary['point_estimate_reference_sha256'])
        check(prefix+'reference stopping flag',summary['reference_optimizer_success']==reference['optimizer_success'])
        check(prefix+'unrestricted-reference flag',summary['reference_is_unrestricted']==(summary['point_estimate_reference_file'] in summary['full_fit_files']))
        close(prefix+'unrestricted reference gap',summary['unrestricted_reference_gap_chi2'],min(x['chi2'] for x in full)-reference['chi2'])
        at=model.labels.index('shift_2382');p=np.array(reference['parameters']);r,J,_=model.evaluate(p)
        close(prefix+'reference D',summary['D_m_s'],1000*p[at])
        information=curvature(J,at)
        close(prefix+'independent local sigma',summary['conditional_sigma_m_s'],information['sigma_qr_m_s'],rtol=1e-10)
        # The refined six-component basin contains a narrow component at a bound.
        # Recheck derivatives at the actual final endpoint, not only old starts.
        for label in ['shift_2382','logb_1','logfwhm_2382']:
            j=model.labels.index(label);step=2e-5 if label.startswith('shift_') else 2e-6
            plus=p.copy();minus=p.copy();plus[j]+=step;minus[j]-=step
            numeric=(model.fun(plus)-model.fun(minus))/(2*step)
            relative=float(np.linalg.norm(numeric-J[:,j])/max(np.linalg.norm(numeric),1e-12))
            check(prefix+'final-basin finite difference '+label,relative<3e-5,relative)
        recorded_set={file for file,obj in records.items() if obj['target']==target}
        new_full_set=set(summary['full_fit_files'])&recorded_set
        new_null_set=set(summary.get('null_fit_files',[]))&recorded_set
        check(prefix+'complete fit inventory',recorded_set==set(summary['all_profile_fit_files'])|new_full_set|new_null_set)
        points=summary['profile_points']
        coordinates=sorted({records[x['fit_file']]['fixed_shift_km_s'] for x in points})
        check(prefix+'unique final fixed-coordinate inventory',len(points)==len(coordinates))
        for point,coordinate in zip(points,coordinates):
            candidates=[x for x in fixed if x['fixed_shift_km_s']==coordinate]
            selected=records[point['fit_file']]
            check(prefix+'at least two starts '+str(coordinate),len(candidates)>=2)
            close(prefix+'minimum of both starts '+str(coordinate),point['chi2'],min(x['chi2'] for x in candidates),rtol=0,atol=1e-11)
            close(prefix+'fixed D '+str(coordinate),point['D_m_s'],1000*coordinate,rtol=0,atol=0)
            close(prefix+'conditional objective rise '+str(coordinate),point['delta_chi2'],selected['chi2']-reference['chi2'])
            check(prefix+'selected stopping flag '+str(coordinate),point['optimizer_success']==selected['optimizer_success'])
            check(prefix+'selected nuisance bounds '+str(coordinate),point['active_bounds']==selected['active_bounds'])
        close(prefix+'plotted grid D',summary['profile_grid_D_m_s'],[x['D_m_s'] for x in points])
        close(prefix+'plotted grid objective rise',summary['profile_delta_chi2'],[x['chi2']-reference['chi2'] for x in points])
        check(prefix+'all selected stopping flag',summary['profile_all_selected_terminated']==all(x['optimizer_success'] for x in points))
        check(prefix+'profile includes D0',any(x['D_m_s']==0 for x in points))
        design_path=folder/('refined_design.json' if refined and (folder/'refined_design.json').exists() else 'design.json')
        design=json.loads(design_path.read_text())
        configured_grid=design.get('grid_km_s',design.get('initial_grid_km_s'))
        check(prefix+'profile includes configured grid',all(any(abs(x['D_m_s']/1000-v)<1e-12 for x in points) for v in configured_grid))
        center=summary['D_m_s'];baseq=reference['chi2']
        for contour in summary['profile_intervals']:
            level=contour['delta_chi2_level']
            check(prefix+'contour closed flag '+str(level),contour['closed']==all(contour[side] is not None for side in ['lower','upper']))
            for side,direction in [('lower',-1),('upper',1)]:
                candidates=sorted([x for x in points if direction*(x['D_m_s']-center)>1e-7],key=lambda x:abs(x['D_m_s']-center))
                inside={'D_m_s':center,'chi2':baseq,'optimizer_success':reference['optimizer_success'],'fit_file':summary['point_estimate_reference_file']}
                outside=None
                for candidate in candidates:
                    if candidate['chi2']-baseq>=level:
                        outside=candidate;break
                    inside=candidate
                given=contour[side]
                check(prefix+f'contour {level} {side} bracket exists', (given is not None)==(outside is not None))
                if outside is None or given is None:continue
                fraction=(baseq+level-inside['chi2'])/(outside['chi2']-inside['chi2'])
                crossing=inside['D_m_s']+fraction*(outside['D_m_s']-inside['D_m_s'])
                close(prefix+f'contour {level} {side} interpolation',given['crossing_m_s'],crossing)
                close(prefix+f'contour {level} {side} bracket',given['bracket_m_s'],sorted([inside['D_m_s'],outside['D_m_s']]))
                check(prefix+f'contour {level} {side} sources',given['fit_files']==[inside['fit_file'],outside['fit_file']])
                check(prefix+f'contour {level} {side} stopped',given['both_endpoints_terminated']==(inside['optimizer_success'] and outside['optimizer_success']))
        if refined:
            check(prefix+'first pass products unchanged',all(sha(ROOT/file)==digest for file,digest in summary['first_pass_sha256'].items()))
            comparison=summary['full_comparison']
            null=json.loads((ROOT/comparison['null_file']).read_text());alternative=json.loads((ROOT/comparison['alternative_file']).read_text())
            nulls=[json.loads((ROOT/file).read_text()) for file in summary['null_fit_files']]
            close(prefix+'lowest conventional endpoint',null['chi2'],min(x['chi2'] for x in nulls))
            close(prefix+'conventional comparison chi2',comparison['chi2_null'],null['chi2'])
            close(prefix+'alternative comparison chi2',comparison['chi2_alternative'],alternative['chi2'])
            close(prefix+'conventional comparison improvement',comparison['delta_chi2'],null['chi2']-alternative['chi2'])
            nullmodel=expected.ExpectedModel(target,summary['ncomp'],False,33)
            check(prefix+'conventional added dimension',comparison['added_parameters']==model.npar-nullmodel.npar)
            lookup=dict(zip(null['labels'],null['parameters']))
            embed=np.array([lookup.get(label,0.) for label in model.labels])
            rn,_,_=nullmodel.evaluate(np.array(null['parameters']))
            ra,_,_=model.evaluate(embed)
            close(prefix+'all-shift-zero H0 nests exactly',rn,ra,rtol=0,atol=0)
            check(prefix+'resolved reference tolerance',summary['profile_baseline_resolved_to_tolerance']==(summary['unrestricted_reference_gap_chi2']<=summary['baseline_tolerance_chi2']))
        rows.append(dict(target=target,D_m_s=summary['D_m_s'],sigma_local_m_s=summary['conditional_sigma_m_s'],reference_is_unrestricted=summary['reference_is_unrestricted'],all_selected_terminated=summary['profile_all_selected_terminated'],n_profile_points=len(points),n_new_fit_products=len(recorded_set),contours=summary['profile_intervals'],full_comparison=summary.get('full_comparison')))
    return rows


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--partial',action='store_true');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    rows=baseline_review()
    probes=actual_wrapper_probes()
    records=fit_products_review()
    profiles=summary_review(records,args.partial)
    aggregate=ROOT/'results/common_pair/refined_summary.json'
    if aggregate.exists():
        top=json.loads(aggregate.read_text())
        check('aggregate contains both targets',len(top['comparisons'])==len(pilot.CONFIG))
        for row in top['comparisons']:
            saved=json.loads((aggregate.parent/row['target']/'refined_summary.json').read_text())
            check(row['target']+' aggregate exact per-target summary',row==saved)
    result=dict(scope='Independent common-pair mathematical/source replay; no optimization or physical significance certification.',
                passed=all(c['passed'] for c in CHECKS),n_checks=len(CHECKS),
                partial=args.partial,n_new_fit_products=len(records),wrapper_probes=probes,
                validator_sha256=sha(__file__),baseline_results=rows,profile_results=profiles,
                failures=[c for c in CHECKS if not c['passed']],checks=CHECKS)
    (OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))


if __name__=='__main__':
    main()
