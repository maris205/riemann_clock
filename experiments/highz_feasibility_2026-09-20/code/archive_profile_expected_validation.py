#!/usr/bin/env python3
"""Independent raw-FITS error-row and frozen-architecture sensitivity audit."""
import hashlib
import io
import json
import sys
import tarfile
from pathlib import Path
import numpy as np

BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'.deps'))
from astropy.io import fits
import archive_profile_pilot as pilot
import archive_profile_expected_noise as sensitivity

OUT=sensitivity.OUT
CHECKS=[]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def check(name,ok,detail=None):CHECKS.append(dict(name=name,passed=bool(ok),detail=detail))
def close(name,a,b,rtol=2e-11,atol=2e-11):
    a,b=np.asarray(a),np.asarray(b)
    ok=a.shape==b.shape and np.allclose(a,b,rtol=rtol,atol=atol,equal_nan=True)
    detail=float(np.nanmax(np.abs(a.astype(float)-b.astype(float)))) if a.shape==b.shape and a.size else None
    check(name,ok,detail)

def main():
    summary=json.loads((OUT/'summary.json').read_text())
    raw={}
    for target in pilot.CONFIG:
        meta=json.loads((pilot.DATA/f'{target}_metadata.json').read_text())
        path=BASE/meta['archive_file']
        check(target+' archive hash',sha(path)==meta['archive_sha256'])
        with tarfile.open(path) as archive:data=archive.extractfile(f'{target}/{target}.fits').read()
        check(target+' raw FITS hash',hashlib.sha256(data).hexdigest()==meta['fits_sha256'])
        with fits.open(io.BytesIO(data)) as hdus:
            raw[target]=hdus[0].data.copy()
            header=hdus[0].header
            check(target+' row2 label',header['ARRAY3']=='Norm. expected fluctuation')
        coadd=np.load(pilot.DATA/f'{target}.npz')
        close(target+' full row2 raw array',coadd['expected_fluctuation'],raw[target][2],rtol=0,atol=0)

    fit_paths=sorted(p for p in OUT.glob('*.json') if '_n' in p.stem and ('_null_' in p.stem or '_shift_' in p.stem))
    fit_records={}
    for path in fit_paths:
        obj=json.loads(path.read_text());fit_records[obj['name']]=obj
        target=obj['target'];n=obj['ncomp'];free=obj['free_shifts'];p=np.asarray(obj['parameters']);prefix=obj['name']
        base=pilot.Model(target,n,free,obj['oversample'])
        model=sensitivity.ExpectedModel(target,n,free,obj['oversample'])
        base_residual,base_jacobian,_=base.evaluate(p)
        residual,jacobian,profiles=model.evaluate(p)
        scaled_r=[];scaled_j=[]
        for ordinary,line,sl in zip(base.lines,model.lines,model.slices):
            key=line['key'];ids=line['source_indices'];good=line['good'];q=prefix+' FeII'+str(key)
            close(q+' same native indices',ids,ordinary['source_indices'],rtol=0,atol=0)
            close(q+' same mask',good,ordinary['good'],rtol=0,atol=0)
            close(q+' source raw row1 statistical error',line['statistical_error'],raw[target][1,ids],rtol=0,atol=0)
            close(q+' source raw row2 expected fluctuation',line['error'],raw[target][2,ids],rtol=0,atol=0)
            check(q+' frozen retained expected errors finite positive',np.all(np.isfinite(line['error'][good])&(line['error'][good]>0)))
            close(q+' raw source flux',line['flux'],raw[target][0,ids],rtol=0,atol=0)
            scale=raw[target][1,ids][good]/raw[target][2,ids][good]
            scaled_r.extend(base_residual[sl]*scale)
            scaled_j.extend(base_jacobian[sl]*scale[:,None])
        close(prefix+' independent residual reweighting',residual,scaled_r)
        close(prefix+' independent analytic Jacobian reweighting',jacobian,scaled_j)
        check(prefix+' frozen labels and architecture',model.labels==base.labels==obj['labels'])
        check(prefix+' same dimension',model.ndata==base.ndata==obj['ndata'] and model.npar==base.npar==obj['npar'])
        check(prefix+' sampling33',obj['oversample']==33)
        check(prefix+' error kind disclosed','row2' in obj['error_kind'])
        check(prefix+' wrapper hash',obj['sensitivity_script_sha256']==sha(sensitivity.__file__))
        check(prefix+' base model hash',obj['base_model_script_sha256']==obj['model_script_sha256']==sha(pilot.__file__))
        check(prefix+' coadd hash',obj['coadd_sha256']==sha(pilot.DATA/f'{target}.npz'))
        check(prefix+' atomic hash',obj['atomic_sha256']==sha(pilot.ATOMIC))
        check(prefix+' source windows hash',all(sha(k)==v for k,v in obj['source_files'].items()))
        check(prefix+' explicit frozen-selection metadata','no component-count reselection' in obj['selection_rule'])
        check(prefix+' termination status consistency',obj['optimizer_success']==(obj['optimizer_status']>0))
        check(prefix+' finite recorded optimality',np.isfinite(obj['optimality']))
        check(prefix+' feasible optimized parameters',np.all(p>=obj['bounds_lower']) and np.all(p<=obj['bounds_upper']))
        close(prefix+' objective',obj['chi2'],residual@residual)
        close(prefix+' objective per nominal ndf',obj['chi2_per_ndf'],(residual@residual)/(model.ndata-model.npar))
        saved=np.load(path.with_suffix('.npz'))
        close(prefix+' saved residual',saved['residual'],residual)
        close(prefix+' saved Jacobian',saved['jacobian'],jacobian)
        close(prefix+' saved parameter vector',saved['parameters'],p,rtol=0,atol=0)
        for line,profile,sl,stated in zip(model.lines,profiles,model.slices,obj['per_line']):
            q=prefix+' FeII'+str(line['key'])
            close(q+' saved model',saved[f'{line["key"]}_model'],profile)
            close(q+' recorded line chi2',stated['chi2'],residual[sl]@residual[sl])
            for k in ['wave','v','flux','error','good','source_indices']:
                close(q+' saved '+k,saved[f'{line["key"]}_{k}'],line[k],rtol=0,atol=0)

    results=[]
    check('two target summaries',len(summary['comparisons'])==2)
    for row in summary['comparisons']:
        target=row['target'];original=json.loads((sensitivity.BASEOUT/f'{target}_selected_refined.json').read_text())
        h0=fit_records[row['expected_null']];h1=fit_records[row['expected_alternative']]
        check(target+' exact primary selection provenance',row['row1_null']==original['null'] and row['row1_alternative']==original['alternative'])
        check(target+' frozen component count',row['ncomp']==original['ncomp']==h0['ncomp']==h1['ncomp'])
        check(target+' frozen data dimensions',row['ndata']==original['ndata']==h0['ndata']==h1['ndata'])
        check(target+' both selected fits terminated successfully',row['both_optimizer_success']==h0['optimizer_success']==h1['optimizer_success']==True)
        check(target+' extra shifts count',row['extra_shifts']==h1['npar']-h0['npar']==original['extra_shifts'])
        close(target+' summary row1 null',row['row1_chi2_null'],original['chi2_null'])
        close(target+' summary row1 alternative',row['row1_chi2_alternative'],original['chi2_alternative'])
        close(target+' summary row1 delta',row['row1_delta_chi2'],original['delta_chi2'])
        close(target+' summary expected null',row['expected_chi2_null'],h0['chi2'])
        close(target+' summary expected alternative',row['expected_chi2_alternative'],h1['chi2'])
        close(target+' summary expected delta',row['expected_delta_chi2'],h0['chi2']-h1['chi2'])
        old0=np.load(sensitivity.BASEOUT/(original['null']+'.npz'));old1=np.load(sensitivity.BASEOUT/(original['alternative']+'.npz'))
        new0=np.load(OUT/(h0['name']+'.npz'));new1=np.load(OUT/(h1['name']+'.npz'))
        for k in original['configuration']['lines']:
            for key in ['source_indices','good','wave','flux']:
                q=f'{k}_{key}'
                close(target+' H0 old to expected '+q,old0[q],new0[q],rtol=0,atol=0)
                close(target+' H1 old to expected '+q,old1[q],new1[q],rtol=0,atol=0)
                close(target+' expected H0 H1 '+q,new0[q],new1[q],rtol=0,atol=0)
        for er in row['per_line_error_ratio']:
            k=er['line'];ids=new0[f'{k}_source_indices'];good=new0[f'{k}_good'];ratio=raw[target][2,ids][good]/raw[target][1,ids][good]
            check(target+' line ratio retained count '+str(k),er['n']==int(good.sum()))
            for desc,func in [('median',np.median),('min',np.min),('max',np.max)]:
                close(target+' line ratio '+desc+' '+str(k),er[desc+'_expected_over_statistical'],func(ratio))
        results.append(dict(target=target,ncomp=row['ncomp'],ndata=row['ndata'],extra_shifts=row['extra_shifts'],row1_delta_chi2=row['row1_delta_chi2'],row2_delta_chi2=row['expected_delta_chi2'],row2_chi2_null=row['expected_chi2_null'],row2_chi2_alternative=row['expected_chi2_alternative'],both_optimizer_success=row['both_optimizer_success']))
    result=dict(scope='Independent original FITS row2, unchanged native masks/architecture, weight/Jacobian/objective replay and selected result bookkeeping; no physical significance certification.',passed=all(c['passed'] for c in CHECKS),n_checks=len(CHECKS),n_saved_fits=len(fit_paths),sensitivity_script_sha256=sha(sensitivity.__file__),validator_sha256=sha(__file__),summary_sha256=sha(OUT/'summary.json'),results=results,failures=[c for c in CHECKS if not c['passed']],checks=CHECKS)
    (OUT/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# Independent expected-fluctuation sensitivity validation','',f'Passed: **{result["passed"]}**; {len(CHECKS)} checks; {len(fit_paths)} fit products replayed.','', 'The error array was reconstructed directly from original SQUAD FITS primary-array row2, independently of the processed NPZ and sensitivity wrapper. Original row1 statistical errors were also checked. Selected native indices, wavelengths, flux, masks, component counts and parameter labels match the earlier row1 analysis. Residuals and analytic Jacobians are independently reproduced by multiplying the original model outputs by raw row1/raw row2 ratios. All source hashes, saved arrays, objectives, summary arithmetic and selected optimizer termination flags were checked.','', '|Target|Components|Pixels|Row1 Δχ²|Row2 H0 χ²|Row2 H1 χ²|Row2 Δχ²|Extra shifts|','|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in results:lines.append(f'|{r["target"]}|{r["ncomp"]}|{r["ndata"]}|{r["row1_delta_chi2"]:.6f}|{r["row2_chi2_null"]:.6f}|{r["row2_chi2_alternative"]:.6f}|{r["row2_delta_chi2"]:.6f}|{r["extra_shifts"]}|')
    lines+=['','This control changes diagonal pixel weights only. It does not model covariance, reselect gas architecture, calibrate wavelength distortions, resolve all blends or establish a global optimum. Recorded successful stopping is a numerical fact, not proof of physical adequacy. Neither relative shifts nor their χ² improvement directly measure alpha or a cosmological time law.']
    if result['failures']:lines+=['','Failed checks:','']+[f'- {c["name"]}: {c["detail"]}' for c in result['failures']]
    (BASE/'reports/archive_profile_expected_validation.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2),flush=True)

if __name__=='__main__':main()
