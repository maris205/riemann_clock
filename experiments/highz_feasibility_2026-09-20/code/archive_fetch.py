#!/usr/bin/env python3
"""Download the five frozen SQUAD archive candidates and decode actual arrays."""
from __future__ import annotations
import concurrent.futures as cf
import hashlib, io, json, sys, tarfile
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import requests
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'.deps'))
from astropy.io import fits
RAW=ROOT/'data/raw/archive_expansion'; DATA=ROOT/'data/processed/archive_expansion'
OUT=ROOT/'results/archive_expansion'
def sha(x): return hashlib.sha256(x).hexdigest()
def plain(x):
    if isinstance(x,bytes): return x.decode('ascii',errors='replace').strip()
    if isinstance(x,np.ndarray): return [plain(i) for i in x]
    if isinstance(x,np.generic): return x.item()
    return x
def fetch_one(selected,resource):
    name=selected['target']; path=RAW/(name+'_Final_Spectrum.tar.gz')
    url=resource['url'].replace('data-portal.hpc.swin.edu.au','dmc.datacentral.org.au')
    if not path.exists():
        response=requests.get(url,timeout=(20,180)); response.raise_for_status()
        if not 1000<len(response.content)<50_000_000: raise ValueError('Unexpected archive size')
        path.write_bytes(response.content)
    archive=path.read_bytes()
    with tarfile.open(fileobj=io.BytesIO(archive),mode='r:gz') as tar:
        item=tar.getmember(f'{name}/{name}.fits')
        if not item.isfile() or item.size>50_000_000: raise ValueError('Unexpected FITS member')
        fb=tar.extractfile(item).read()
        upl=[m for m in tar.getmembers() if m.isfile() and m.name.endswith('.upl')]
        if upl: (RAW/(name+'.upl')).write_bytes(tar.extractfile(upl[0]).read())
    with fits.open(io.BytesIO(fb),memmap=False) as hdul:
        a=np.asarray(hdul[0].data,dtype=np.float64); h=hdul[0].header
        assert a.ndim==2 and a.shape[0]==9 and h['DC-FLAG']==1
        wl=10**(h['CRVAL1']+(np.arange(a.shape[1])+1-h['CRPIX1'])*h['CD1_1'])
        assert np.all(np.diff(wl)>0)
        assert np.isclose(wl[0],h['UP_WLSRT'],rtol=1e-10)
        assert np.isclose(wl[-1],h['UP_WLEND'],rtol=1e-10)
        status=np.rint(a[4]).astype(np.int16)
        valid=(status==1)&np.isfinite(a[0])&np.isfinite(a[1])&(a[1]>0)
        np.savez_compressed(DATA/(name+'.npz'),wavelength_AA=wl,flux=a[0],error=a[1],
            expected_fluctuation=a[2],continuum=a[3],status=status,valid=valid,
            n_contrib_before=a[5],n_contrib_after=a[6],chi2_before=a[7],chi2_after=a[8])
        tables={str(i):{'name':ext.name,'columns':list(ext.columns.names),
            'rows':[{col:plain(row[col]) for col in ext.columns.names} for row in ext.data]}
            for i,ext in enumerate(hdul[1:],1)}
        dates=[r['UTDate'] for r in tables['3']['rows']]
        meta={**selected,'archive_url':url,'archive_file':str(path.relative_to(ROOT)),
            'archive_bytes':len(archive),'archive_sha256':sha(archive),'fits_sha256':sha(fb),
            'access_date_utc':datetime.now(timezone.utc).isoformat(),
            'pixels':a.shape[1],'valid_pixels':int(valid.sum()),'native_dispersion_km_s':float(h['UP_DISP']),
            'n_exposures':int(h['UP_NEXP']),'exposure_s':float(h['UP_TEXP']),
            'observation_date_range':[min(dates),max(dates)],
            'wavelength_AA_range':[float(wl[0]),float(wl[-1])],
            'header':dict(h),'tables':tables}
    (DATA/(name+'_metadata.json')).write_text(json.dumps(meta,indent=2,default=str)+'\n')
    return {k:v for k,v in meta.items() if k not in ['header','tables']}
def main():
    RAW.mkdir(parents=True,exist_ok=True); DATA.mkdir(parents=True,exist_ok=True)
    screen=json.loads((OUT/'catalogue_screen.json').read_text())
    resources=json.loads((ROOT/'data/raw/catalog_portal_metadata.json').read_text())['result']['resources']
    lookup={r['name']:r for r in resources}
    results=[]
    with cf.ThreadPoolExecutor(max_workers=5) as executor:
        future={executor.submit(fetch_one,s,lookup[s['target']+'_Final_Spectrum.tar.gz']):s for s in screen['selected']}
        for f in cf.as_completed(future):
            row=f.result();results.append(row)
            print(json.dumps({k:row[k] for k in ['target','z_abs','archive_bytes','pixels','n_exposures']}),flush=True)
    results.sort(key=lambda r:r['target'])
    (OUT/'source_manifest.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),
        'selection_policy_sha256':sha((OUT/'selection_policy.json').read_bytes()),
        'catalogue_screen_sha256':sha((OUT/'catalogue_screen.json').read_bytes()),
        'sources':results,'license':'Spectrum portal CC BY-SA, version unspecified; metadata/code CC BY 4.0. Cite Murphy et al.2019 DOI10.1093/mnras/sty2834.',
        'processing':'Decoded published continuum-normalized coadds. No rebinning or raw-photon extraction. These are historical archival observations, not new 2026 acquisitions.'},indent=2)+'\n')
if __name__=='__main__': main()
