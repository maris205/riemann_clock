#!/usr/bin/env python3
"""Restore the published spectral snapshot, with pinned revision and SHA-256 checks.

Uses Python's standard library. Public downloads require no Hugging Face token.
Existing files are verified and preserved; mismatches stop rather than overwrite.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import fnmatch
import hashlib
import json
from pathlib import Path, PurePosixPath
import tempfile
import time
from urllib.parse import quote
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / 'data/huggingface_release.json'
EXPERIMENT = 'experiments/highz_feasibility_2026-09-20/data/'


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def public_url(release, path):
    return ('https://huggingface.co/datasets/' + release['repo_id'] + '/resolve/'
            + release['revision'] + '/' + quote(path, safe='/') + '?download=true')


def request(url):
    for attempt in range(4):
        try:
            return urlopen(Request(url, headers={'User-Agent': 'riemann-clock-reproduction/1'}),
                           timeout=180)
        except HTTPError as exc:
            if exc.code not in {408, 429, 500, 502, 503, 504} or attempt == 3:
                raise
        except (URLError, OSError):
            if attempt == 3:
                raise
        time.sleep(1 + attempt)


def load_manifest(release):
    with request(public_url(release, 'manifest.json')) as response:
        content = response.read(10 * 1024 * 1024 + 1)
    if len(content) > 10 * 1024 * 1024:
        raise ValueError('Unexpectedly large dataset manifest')
    if hashlib.sha256(content).hexdigest() != release['manifest_sha256']:
        raise ValueError('Dataset manifest does not match the pinned release')
    return json.loads(content)


def destination(base, row):
    rel = PurePosixPath(row['project_relative_path'])
    if rel.is_absolute() or '..' in rel.parts:
        raise ValueError('Unsafe destination in manifest')
    prefix = EXPERIMENT + ('raw/' if row['category'] == 'source_spectrum' else 'processed/')
    if not str(rel).startswith(prefix):
        raise ValueError('Destination is outside the selected spectral data directories')
    path = base.joinpath(*rel.parts)
    if not path.resolve().is_relative_to(base.resolve()):
        raise ValueError('Destination resolves outside the project directory')
    return path


def check_file(path, row):
    if path.stat().st_size != row['bytes'] or sha256(path) != row['sha256']:
        raise ValueError(f'Existing file differs from the published snapshot: {path}')


def restore_one(base, release, row, verify_only):
    path = destination(base, row)
    if path.exists():
        check_file(path, row)
        return 'already_verified'
    if verify_only:
        raise FileNotFoundError(f'Missing required file: {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    last_error = None
    for attempt in range(3):
        temporary = None
        try:
            h = hashlib.sha256()
            size = 0
            with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.' + path.name + '.',
                                             suffix='.download', delete=False) as output:
                temporary = Path(output.name)
                with request(public_url(release, row['dataset_path'])) as response:
                    while chunk := response.read(4 * 1024 * 1024):
                        size += len(chunk)
                        if size > row['bytes']:
                            raise ValueError('Download exceeds manifest size')
                        h.update(chunk)
                        output.write(chunk)
            if size != row['bytes'] or h.hexdigest() != row['sha256']:
                raise ValueError('Downloaded bytes do not match the pinned manifest')
            # A concurrent invocation must not silently replace an existing file.
            try:
                path.hardlink_to(temporary)
            except FileExistsError:
                check_file(path, row)
                return 'already_verified'
            return 'downloaded_verified'
        except Exception as exc:
            last_error = exc
            if attempt < 2:
                time.sleep(1 + attempt)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    raise RuntimeError(f'Could not restore {row["dataset_path"]}: {last_error}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, default=ROOT,
                        help='Project root to restore into (default: this checkout)')
    parser.add_argument('--include-processed', action='store_true',
                        help='Also restore the processed arrays already included in Git')
    parser.add_argument('--only', help='Select dataset paths matching this shell-style pattern')
    parser.add_argument('--verify-only', action='store_true',
                        help='Verify local files against the public manifest; download no spectra')
    parser.add_argument('--dry-run', action='store_true',
                        help='Fetch and verify the small manifest; list selected files without changes')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.workers <= 16:
        parser.error('--workers must be between 1 and 16')
    release = json.loads(RELEASE.read_text())
    manifest = load_manifest(release)
    categories = {'source_spectrum', 'processed'} if args.include_processed else {'source_spectrum'}
    rows = [r for r in manifest['files'] if r['category'] in categories
            and (not args.only or fnmatch.fnmatchcase(r['dataset_path'], args.only))]
    if not rows:
        parser.error('No spectral files match the selection')
    paths = [destination(args.destination, row) for row in rows]
    if len(paths) != len(set(paths)):
        raise ValueError('Duplicate restore destinations in manifest')
    total = sum(r['bytes'] for r in rows)
    if args.dry_run:
        print(json.dumps({'repo_id': release['repo_id'], 'revision': release['revision'],
                          'selected_files': len(rows), 'selected_bytes': total,
                          'dataset_paths': [r['dataset_path'] for r in rows]}, indent=2))
        return
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        statuses = list(executor.map(lambda r: restore_one(args.destination, release, r, args.verify_only), rows))
    print(json.dumps({'status': 'PASS', 'repo_id': release['repo_id'], 'revision': release['revision'],
                      'verified_files': len(rows), 'verified_bytes': total,
                      'downloaded_verified': statuses.count('downloaded_verified'),
                      'already_verified': statuses.count('already_verified')}, indent=2))


if __name__ == '__main__':
    main()
