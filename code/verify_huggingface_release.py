#!/usr/bin/env python3
"""Verify anonymous access and every published file's identity at the pinned revision.

Large files are checked against the Hub's SHA-256 metadata; ordinary Git files
are downloaded and hashed. This does not rerun the scientific analysis.
"""
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json

from fetch_huggingface_spectra import ROOT, RELEASE, load_manifest, public_url, request


def main():
    release = json.loads(RELEASE.read_text())
    manifest = load_manifest(release)
    url = ('https://huggingface.co/api/datasets/' + release['repo_id']
           + '/revision/' + release['revision'] + '?blobs=true')
    with request(url) as response:
        info = json.load(response)
    checks = []

    def check(name, passed, detail=None):
        row = {'name': name, 'passed': bool(passed)}
        if detail is not None:
            row['detail'] = detail
        checks.append(row)

    check('Anonymous metadata identifies the pinned public dataset',
          info['id'] == release['repo_id'] and info['sha'] == release['revision']
          and info.get('private') is False and not info.get('gated'))
    siblings = {r['rfilename']: r for r in info['siblings']}
    files = manifest['files']
    expected = {r['dataset_path'] for r in files}
    check('Remote file inventory equals the manifest plus manifest and Git attributes',
          expected | {'manifest.json'} <= set(siblings)
          and not (set(siblings) - expected - {'manifest.json', '.gitattributes'}))
    check('Manifest paths are unique', len(expected) == len(files))
    methods = {'server_sha256': 0, 'anonymous_download_sha256': 0}
    def verify_file(row):
        name = row['dataset_path']
        remote = siblings.get(name, {})
        lfs = remote.get('lfs') or {}
        size = remote.get('size', lfs.get('size'))
        if lfs:
            digest = lfs.get('sha256', lfs.get('oid'))
            method = 'server_sha256'
        else:
            with request(public_url(release, name)) as response:
                content = response.read(row['bytes'] + 1)
            digest = hashlib.sha256(content).hexdigest()
            size = len(content)
            method = 'anonymous_download_sha256'
        return name, size == row['bytes'] and digest == row['sha256'], method
    with ThreadPoolExecutor(max_workers=4) as executor:
        for name, passed, method in executor.map(verify_file, files):
            methods[method] += 1
            check('Published file identity: ' + name, passed)
    for category, count, size in [('source_spectrum', 52, 1366293902),
                                  ('processed', 424, 260049012)]:
        rows = [r for r in files if r['category'] == category]
        check('Frozen ' + category + ' inventory',
              len(rows) == count and sum(r['bytes'] for r in rows) == size)
    output = {
        'scope': __doc__.strip(), 'generated_utc': datetime.now(timezone.utc).isoformat(),
        'repo_id': release['repo_id'], 'revision': release['revision'],
        'manifest_sha256': release['manifest_sha256'],
        'passed': all(c['passed'] for c in checks),
        'checks_passed': sum(c['passed'] for c in checks), 'checks_total': len(checks),
        'published_files_checked': len(files), 'verification_methods': methods,
        'checks': checks,
    }
    path = ROOT / 'reports/huggingface_release_verification.json'
    path.write_text(json.dumps(output, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({k: v for k, v in output.items() if k != 'checks'}, indent=2))
    if not output['passed']:
        print(json.dumps([r for r in checks if not r['passed']], indent=2))
        raise SystemExit(1)


if __name__ == '__main__':
    main()
