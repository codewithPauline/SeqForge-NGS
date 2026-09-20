"""Shared verification and interval helpers for benchmark preparation."""
import hashlib
import json
import re
import urllib.request
from pathlib import Path

REQUIRED = {'reference/reference.fa', 'reference/reference.fa.fai', 'reference/reference.dict',
            'reference/targets.bed', 'reads.bam', 'reads.bam.bai', 'truth.vcf.gz', 'truth.vcf.gz.tbi', 'evaluation.bed'}


def checksum(path, algorithm='sha256'):
    digest = hashlib.new(algorithm)
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def validate_manifest(m):
    if m.get('schema_version') != 1 or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}', m.get('sample', '')):
        raise ValueError('Unsupported manifest schema or invalid sample ID')
    region = m['scoring_region']
    if region.get('coordinate_system') != '0-based half-open' or not re.fullmatch(r'[A-Za-z0-9_.]+', region['chrom']):
        raise ValueError('Invalid scoring coordinate system or contig')
    if not all(type(region[k]) is int for k in ('start', 'end')) or not 0 <= region['start'] < region['end']:
        raise ValueError('Invalid scoring coordinates')
    if type(m['padding_bases']) is not int or m['padding_bases'] < 0:
        raise ValueError('Padding must be a nonnegative integer')
    if type(m['min_depth']) is not int or m['min_depth'] < 1 or type(m['min_qual']) not in (int, float) or not 0 <= m['min_qual'] < float('inf'):
        raise ValueError('Invalid QUAL/depth thresholds')


def download(asset, directory):
    name = asset['filename']
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]+', name) or not asset['url'].startswith('https://'):
        raise ValueError('Asset requires a safe basename and HTTPS URL')
    if not re.fullmatch('[0-9a-f]{32}', asset['md5']) or type(asset['bytes']) is not int or asset['bytes'] < 1:
        raise ValueError('Asset requires provider MD5 and positive byte size')
    path = Path(directory) / name
    if path.exists():
        if path.stat().st_size != asset['bytes'] or checksum(path, 'md5') != asset['md5']:
            raise ValueError(f'Cached asset failed verification: {path}; remove it before retrying')
        print(f'Verified cache: {name}', flush=True)
        return path
    print(f'Downloading {name} ({asset["bytes"] / 1e6:.1f} MB)', flush=True)
    temporary = path.with_name(name + '.part')
    digest, size = hashlib.md5(), 0
    try:
        with urllib.request.urlopen(asset['url'], timeout=120) as response, temporary.open('wb') as target:
            while block := response.read(1024 * 1024):
                size += len(block)
                if size > asset['bytes']:
                    raise ValueError(f'Asset exceeds declared size: {name}')
                digest.update(block)
                target.write(block)
        if size != asset['bytes'] or digest.hexdigest() != asset['md5']:
            raise ValueError(f'Asset size or MD5 mismatch: {name}')
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    print(f'Verified download: {name}', flush=True)
    return path


def intersect_bed(path, chrom, start, end):
    intervals = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.split()
        if len(fields) < 3:
            raise ValueError('Invalid confident BED record')
        if fields[0] != chrom:
            continue
        left, right = map(int, fields[1:3])
        if not 0 <= left < right:
            raise ValueError('Invalid confident BED coordinates')
        if max(start, left) < min(end, right):
            intervals.append((max(start, left), min(end, right)))
    merged = []
    for left, right in sorted(intervals):
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(right, merged[-1][1]))
        else:
            merged.append((left, right))
    if not merged:
        raise ValueError('Scoring region has no overlap with confident regions')
    return merged


def inventory(bundle):
    bundle = Path(bundle)
    return {name: {'bytes': (bundle / name).stat().st_size, 'sha256': checksum(bundle / name)} for name in sorted(REQUIRED)}


def verify(bundle):
    bundle = Path(bundle)
    prepared = json.loads((bundle / 'prepared.json').read_text())
    if prepared.get('schema_version') != 1 or set(prepared['files']) != REQUIRED:
        raise ValueError('Unsupported prepared schema or incomplete inventory')
    validate_manifest(prepared['manifest'])
    for name, expected in prepared['files'].items():
        path = bundle / name
        if not path.is_file() or path.stat().st_size != expected['bytes'] or checksum(path) != expected['sha256']:
            raise ValueError(f'Prepared input failed SHA-256 verification: {name}')
    return {'verified': True, 'files': len(REQUIRED), 'sample': prepared['manifest']['sample'], 'prepared_sha256': checksum(bundle / 'prepared.json')}
