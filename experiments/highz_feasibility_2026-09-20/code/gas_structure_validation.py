#!/usr/bin/env python3
"""Independent numerical adapter audit; does not run or alter scientific fits."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import json, hashlib, math, tempfile, datetime, argparse
import numpy as np
from scipy.linalg import solve
from scipy.special import wofz
from espresso_conventional_null import Model, ROOT, C, ZREF, TAU_CONSTANT, AUTHOR_CONT, AUTHOR_ZERO
from espresso_empirical_gls import GLSModel, CONTROL
from gas_structure_refit import MergedGLSModel, merged_components, convert_start

OUT=ROOT/'results/gas_structure'
checks=[]

def check(name, value, detail=None):
    checks.append(dict(name=name,passed=bool(value),detail=detail))
    if not value: print('FAIL',name,detail,flush=True)

def close(name,a,b,atol=1e-9,rtol=1e-10):
    a=np.asarray(a);b=np.asarray(b)
    err=float(np.max(np.abs(a-b))) if a.size else 0.
    check(name,np.allclose(a,b,atol=atol,rtol=rtol),dict(max_absolute_error=err,atol=atol,rtol=rtol))

def independent_groups(comp,threshold):
    order=sorted(range(len(comp)),key=lambda i:comp[i,1]); groups=[]
    for i in order:
        if not groups or comp[i,1]-comp[groups[-1][-1],1]>=threshold: groups.append([i])
        else: groups[-1].append(i)
    return groups

def independent_merge(comp,groups):
    out=[]
    for g in groups:
        total=math.fsum(10.**float(comp[i,0]) for i in g)
        mean=math.fsum(10.**float(comp[i,0])*float(comp[i,1]) for i in g)/total
        var=math.fsum(10.**float(comp[i,0])*(float(comp[i,2])**2/2+(float(comp[i,1])-mean)**2) for i in g)/total
        out.append([math.log10(total),mean,math.sqrt(2*var)])
    return np.array(out)

def expanded_parameters(m,p,full):
    """Replicate each merged cloud as original-group-sized identical subclouds."""
    q=full.initial()[0];lookup=dict(zip(m.labels,p))
    n=full.ncomp
    for j,g in enumerate(m.groups):
        for i in g:
            q[full.labels.index(f'logN_{i}')]=p[j]-math.log10(len(g))
            q[full.labels.index(f'v_{i}')]=p[m.ncomp+j]
            q[full.labels.index(f'logb_{i}')]=p[2*m.ncomp+j]
    for j,k in enumerate(full.labels[full.ngas:],full.ngas):q[j]=lookup[k]
    return q

def fold_jacobian(m,full,J):
    columns=[]
    for family in ['logN','v','logb']:
        for g in m.groups:
            columns.append(np.sum(J[:,[full.labels.index(f'{family}_{i}') for i in g]],axis=1))
    columns += [J[:,full.labels.index(label)] for label in m.labels[m.ngas:]]
    return np.column_stack(columns)

def direct_profile(m,p,L):
    """Independent scalar optical depth and explicit discrete LSF convolution."""
    lookup=dict(zip(m.labels,p)); grid=L['grid'];tau=np.zeros(len(grid))
    line_shift=lookup.get(f'shift_{L["key"]}',0.)
    for i in range(m.ncomp):
        N=10**lookup[f'logN_{i}'];vel=lookup[f'v_{i}'];b=np.exp(lookup[f'logb_{i}'])
        for lam,f,gamma,_,_ in L['atom']:
            u=C*(-np.expm1((vel+line_shift+C*np.log(lam/L['ref'])-grid)/C))/b
            a=gamma*lam*1e-13/(4*np.pi*b)
            tau+=TAU_CONSTANT*N*f*lam/b*wofz(u+1j*a).real
    transmission=np.exp(-tau)
    sigma=L['fwhm']/2.354820045/(m.dv/m.oversample)
    radius=int(6*sigma+.5);x=np.arange(-radius,radius+1)
    kernel=np.exp(-.5*(x/sigma)**2);kernel/=kernel.sum()
    convolved=np.convolve(np.pad(transmission,radius,mode='edge'),kernel,mode='valid')
    sampled=convolved.reshape(-1,m.oversample).mean(axis=1)[L['pad']:-L['pad']]
    k=L['key'];cont=AUTHOR_CONT[k]*(1+lookup.get(f'continuum0_{k}',0)+lookup.get(f'continuum1_{k}',0)*L['x'])
    return cont*(AUTHOR_ZERO[k]+(1-AUTHOR_ZERO[k])*sampled)

def adapter_tests():
    m=MergedGLSModel(free_shifts=True);n=MergedGLSModel(free_shifts=False)
    groups=independent_groups(m.published_comp,m.threshold)
    check('fixed groups reconstruct independently',groups==m.groups)
    check('40 components at prescribed threshold',m.ncomp==40)
    check('each original component assigned exactly once',sorted(sum(groups,[]))==list(range(45)))
    close('independent moment-conserving initializer',m.comp,independent_merge(m.published_comp,groups),2e-13)
    close('total column conserved',np.sum(10**m.comp[:,0]),np.sum(10**m.published_comp[:,0]),rtol=3e-15,atol=0)
    for j,g in enumerate(groups):
        old=m.published_comp[g];N=10**old[:,0];w=N/N.sum();v=m.comp[j,1]
        close(f'group {j} centroid',v,w@old[:,1],1e-12)
        close(f'group {j} Doppler core variance',m.comp[j,2]**2/2,w@(old[:,2]**2/2+(old[:,1]-v)**2),2e-11)
    p,lo,hi=m.initial()
    close('velocity lower bounds stay at published merged centers',lo[m.ncomp:2*m.ncomp],m.comp[:,1]-2)
    close('velocity upper bounds stay at published merged centers',hi[m.ncomp:2*m.ncomp],m.comp[:,1]+2)
    check('null/alternative dimension difference five',m.npar-n.npar==5)
    check('no duplicate parameter labels',len(set(m.labels))==m.npar)
    for L,L0,ch,ch0 in zip(m.lines,n.lines,m.cholesky,n.cholesky):
        k=L['key'];check(f'{k} same null/alternative mask',np.array_equal(L['good'],L0['good']))
        close(f'{k} same null/alternative error',L['error'],L0['error'])
        close(f'{k} same null/alternative GLS whitening',ch,ch0)
    identity=MergedGLSModel(threshold=0,free_shifts=True);original=GLSModel(free_shifts=True)
    qi=expanded_parameters(identity,identity.initial()[0],original)
    ri,Ji,Pi=identity.evaluate(identity.initial()[0]);ro,Jo,Po=original.evaluate(qi)
    close('zero threshold has original 45 component residual',ri,ro,2e-10)
    close('zero threshold has original Jacobian after permutation',Ji,fold_jacobian(identity,original,Jo),3e-9)
    full=GLSModel(free_shifts=True)
    rng=np.random.default_rng(2191040)
    for iteration in range(2):
        trial=p.copy()
        if iteration:
            trial[:m.ncomp]+=rng.normal(0,.005,m.ncomp)
            trial[m.ncomp:2*m.ncomp]+=rng.normal(0,.025,m.ncomp)
            trial[2*m.ncomp:3*m.ncomp]+=rng.normal(0,.005,m.ncomp)
            trial[m.ngas:m.ngas+m.ncont]+=rng.normal(0,.001,m.ncont)
            trial[m.shift_offset:]+=rng.normal(0,.025,len(m.shift_keys))
        r,J,P=m.evaluate(trial);q=expanded_parameters(m,trial,full);rf,Jf,Pf=full.evaluate(q)
        close(f'expansion {iteration} exact residual',r,rf,5e-10)
        close(f'expansion {iteration} exact folded Jacobian',J,fold_jacobian(m,full,Jf),1e-8)
        for L,a,b in zip(m.lines,P,Pf):close(f'expansion {iteration} line {L["key"]} flux',a,b,2e-12)
        if iteration:
            for L,flux in zip(m.lines,P):
                close(f'independent scalar opacity and discrete LSF {L["key"]}',flux,direct_profile(m,trial,L),atol=3e-10,rtol=2e-10)
            # Reconstruct C on real detector indices; compare quadratic and normal matrix.
            rd,Jd,_=Model.evaluate(m,trial)
            for L,sl,ch in zip(m.lines,m.slices,m.cholesky):
                idx=L['idx'][L['good']];lag=abs(np.subtract.outer(idx,idx))
                acf=np.array(json.loads(CONTROL.read_text())['primary']['acf'])
                weight=np.arange(len(acf),0,-1)/len(acf);kernel=acf*weight
                Cij=np.where(lag<len(kernel),kernel[np.minimum(lag,len(kernel)-1)],0.)
                close(f'{L["key"]} covariance exact pixel separation',ch@ch.T,Cij,3e-15)
                direct=solve(Cij,rd[sl],assume_a='pos')
                close(f'{L["key"]} GLS quadratic direct inverse',r[sl]@r[sl],rd[sl]@direct,atol=2e-8,rtol=1e-12)
                # Representative columns include all nuisance families and shift columns.
                cols=[0,m.ncomp,2*m.ncomp,m.ngas,m.ngas+1,m.npar-1]
                J0=Jd[sl][:,cols];Jw=J[sl][:,cols]
                close(f'{L["key"]} GLS Jacobian normal matrix',Jw.T@Jw,J0.T@solve(Cij,J0,assume_a='pos'),atol=2e-7,rtol=1e-11)
            groups_nontrivial=[i for i,g in enumerate(groups) if len(g)>1]
            gas_indices=sorted(set([0,1,20,m.ncomp-1]+groups_nontrivial))
            cols=[j+offset for offset in (0,m.ncomp,2*m.ncomp) for j in gas_indices]+list(range(m.ngas,m.npar))
            for j in cols:
                step=2e-6 if m.labels[j].startswith('v_') or m.labels[j].startswith('shift_') else 1e-6
                high=trial.copy();low=trial.copy();high[j]+=step;low[j]-=step
                fd=(m.evaluate(high)[0]-m.evaluate(low)[0])/(2*step)
                rel=np.linalg.norm(fd-J[:,j])/max(np.linalg.norm(J[:,j]),1e-8)
                check(f'finite difference {m.labels[j]}',rel<3e-4,dict(relative_L2_error=float(rel),step=step))
    # Mapping from an adversarially reordered fitted centroid state must retain fixed membership.
    fullp=original.initial()[0];fullp[45:90]=fullp[45:90][::-1]
    old=dict(labels=original.labels,parameters=fullp.tolist())
    with tempfile.TemporaryDirectory() as td:
        path=Path(td)/'start.json';path.write_text(json.dumps(old));mapped,_,_=convert_start(m,path)
    oldcomp=np.array([[fullp[i],fullp[45+i],np.exp(fullp[90+i])] for i in range(45)])
    expected=independent_merge(oldcomp,groups)
    expectedp=p.copy();expectedp[:40]=expected[:,0];expectedp[40:80]=expected[:,1];expectedp[80:120]=np.log(expected[:,2]);expectedp=np.clip(expectedp,lo+1e-9,hi-1e-9)
    close('45 component mapping uses fixed published groups even if fitted order changes',mapped,expectedp,1e-12)
    check('adversarial fitted centroids actually imply different dynamic groups',independent_groups(oldcomp,m.threshold)!=groups)
    check('source fit initializers never mutate published component array',np.array_equal(m.published_comp,original.comp))
    return m

def saved_result_tests():
    outputs=[]
    for path in sorted(OUT.glob('merged40_*.json')):
        if path.stem.endswith('_checkpoint'):continue
        d=json.loads(path.read_text())
        if 'optimizer_success' not in d:continue
        m=MergedGLSModel(threshold=d['threshold_km_s'],free_shifts=d['free_shifts']);p=np.array(d['parameters']);r,J,P=m.evaluate(p)
        label=path.stem;arrays=np.load(path.with_suffix('.npz'))
        close(f'{label} stored chi2 replay',d['chi2'],r@r,atol=1e-8)
        close(f'{label} saved residual replay',arrays['residuals'],r,atol=1e-10)
        close(f'{label} saved Jacobian replay',arrays['jacobian'],J,atol=1e-8)
        check(f'{label} exact labels',d['labels']==m.labels)
        check(f'{label} exact groups',d['merge_groups']==m.groups)
        p0,lo,hi=m.initial();check(f'{label} final parameters obey published-centered bounds',np.all(p>=lo-1e-12) and np.all(p<=hi+1e-12))
        start=ROOT/d['start_file'] if d['start_file'] else None
        initial,_,_=convert_start(m,start)
        close(f'{label} start replay fixed membership',d['initial_parameters'],initial,atol=1e-11)
        check(f'{label} data dimensions',d['ndata']==m.ndata and d['npar']==m.npar)
        check(f'{label} core code hash',d['core_source_sha256']==hashlib.sha256((ROOT/'code/espresso_conventional_null.py').read_bytes()).hexdigest())
        check(f'{label} fit code hash',d['code_sha256']==hashlib.sha256((ROOT/'code/gas_structure_refit.py').read_bytes()).hexdigest())
        check(f'{label} covariance source hash',d['covariance_source_sha256']==hashlib.sha256(CONTROL.read_bytes()).hexdigest())
        check(f'{label} optimizer success agrees stopping message',d['optimizer_success']==('maximum number' not in d['message'].lower()))
        close(f'{label} covariance kernel',d['covariance_kernel'],m.kernel,atol=1e-15)
        for L,prof in zip(m.lines,P):
            k=L['key'];close(f'{label} {k} profile',arrays[f'{k}_model'],prof,atol=1e-12)
            check(f'{label} {k} mask',np.array_equal(arrays[f'{k}_good'],L['good']))
            close(f'{label} {k} errors',arrays[f'{k}_error'],L['error'],atol=0,rtol=0)
        outputs.append({k:d[k] for k in ['name','chi2','optimizer_success','message','nfev','optimality']})
    comparison=OUT/'comparison.json'
    if comparison.exists():
        s=json.loads(comparison.read_text());a=json.loads((OUT/(s['null']+'.json')).read_text());b=json.loads((OUT/(s['alternative']+'.json')).read_text())
        close('selected delta chi2 arithmetic',s['delta_chi2'],a['chi2']-b['chi2'])
        attempts=[json.loads((OUT/(item['name']+'.json')).read_text()) for item in s['all_final_attempts']]
        for free,chosen in [(False,a),(True,b)]:
            valid=[item for item in attempts if item['free_shifts']==free and item['optimizer_success']]
            check(f'selected {"alternative" if free else "null"} lowest converged retained endpoint',bool(valid) and chosen['chi2']==min(item['chi2'] for item in valid))
        for item,meta in zip(attempts,s['all_final_attempts']):
            check(f'selected attempt {item["name"]} metadata',all(item[key]==meta[key] for key in meta))
        check('selected pair identical ndata',a['ndata']==b['ndata']==s['ndata'])
        check('selected pair exactly five added parameters',b['npar']-a['npar']==s['added_parameters']==5)
        check('selected convergence statuses truthful',s['null_success']==a['optimizer_success'] and s['alternative_success']==b['optimizer_success'])
        close('selected pair equal covariance kernels',a['covariance_kernel'],b['covariance_kernel'])
        for key in s['shifts_m_s']:close(f'selected shift {key}',s['shifts_m_s'][key],b['shifts_m_s'][key])
    return outputs

def source_fingerprints():
    files=[ROOT/'code/gas_structure_refit.py',ROOT/'code/espresso_conventional_null.py',ROOT/'code/espresso_empirical_gls.py',CONTROL,ROOT/'results/gas_structure/design.json']
    files += [ROOT/'data/raw/espresso_null'/name for name in ['fit_r_iso.f18','fit_r_iso.f13','MM_VPFIT_2013-11-10.dat','hes0515m4414.fits']]
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--saved-only',action='store_true');args=parser.parse_args()
    fingerprints=source_fingerprints()
    if args.saved_only:
        previous=json.loads((OUT/'independent_validation.json').read_text())
        if previous.get('source_fingerprints')!=fingerprints:
            raise RuntimeError('Adapter/data inputs changed or lack cached fingerprints; run full audit.')
        checks.extend(c for c in previous['checks'] if not c['name'].startswith(('merged40_','selected ')))
    else:adapter_tests()
    adapter_check_count=len(checks)
    fit_results=saved_result_tests()
    record=dict(source_fingerprints=fingerprints,adapter_check_count=adapter_check_count,adapter_checks_reused=bool(args.saved_only),generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),passed=all(c['passed'] for c in checks),checks_total=len(checks),checks_passed=sum(c['passed'] for c in checks),fits_replayed=fit_results,comparison_present=(OUT/'comparison.json').exists(),checks=checks,scope='Independent numerical checks of fixed group adapter, bounds, initialization, flux/Jacobian equivalence, direct profile construction, finite differences, covariance, and saved fit replay. A pass does not certify global optimizer convergence or scientific adequacy.')
    (OUT/'independent_validation.json').write_text(json.dumps(record,indent=2)+'\n')
    lines=['# Independent audit of merged gas-structure sensitivity','',f'Generated UTC: {record["generated_utc"]}', '',f'Numerical result: **{record["checks_passed"]}/{record["checks_total"]} PASS**.', '', 'The adapter uses fixed groups of the published 45 centroids. The 1.025 km/s adjacency rule produces 40 components. Independent scalar sums reproduce total column, column-weighted velocity and the finite Doppler-core second moment; no finite Lorentzian/Voigt moment is asserted. Bounds remain centered on this published merged initializer.', '', 'Tests include exact expansion of every merged cloud into coincident subclouds in the unchanged 45-component engine, equality of profiles and folded analytic Jacobians, the no-merging limit, scalar optical-depth and explicit Gaussian convolution calculations, finite differences of gas/continuum/shift columns, and direct covariance inversion using true pixel separations. An adversarial input reverses fitted centroid order to verify that initialization does not regroup fitted centroids.', '', f'{len(fit_results)} saved optimizer endpoints replayed. Final comparison present: {record["comparison_present"]}. Numerical replay is distinct from optimizer convergence. All retained statuses and stopping messages are listed below.', '', '| Endpoint | chi2 | Success | nfev | Optimality |','|---|---:|---|---:|---:|']
    for d in fit_results:lines.append(f'| {d["name"]} | {d["chi2"]:.8f} | {d["optimizer_success"]} | {d["nfev"]} | {d["optimality"]:.6g} |')
    lines += ['', 'Scientific limits: this is an exploratory sensitivity protocol frozen before its new refits, not a prospectively preregistered discovery test. One predetermined component-merging architecture cannot exhaust gas-structure uncertainty. The changed architecture also changes centroid bounds through its initializer. A formally successful local least-squares termination is not proof of a global minimum. Empirical continuum covariance transfer and wavelength calibration remain conditional. Comparisons with the 45-component model are a combined structure/initialization sensitivity, not a nested significance test or evidence for a cosmic-time law.', '', 'Full checks: `results/gas_structure/independent_validation.json`.']
    (ROOT/'reports/gas_structure_independent_review.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='checks'},indent=2))
    if not record['passed']:raise SystemExit(1)

if __name__=='__main__':main()
