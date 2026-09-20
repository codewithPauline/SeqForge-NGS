#!/usr/bin/env python3
"""Download pinned HG002 resources and prepare a padded regional calling benchmark."""
import argparse
import concurrent.futures
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from benchmark_io import download, intersect_bed, inventory, validate_manifest


def prepare(manifest_path, output):
    manifest = json.loads(Path(manifest_path).read_text())
    validate_manifest(manifest)
    output = Path(output).resolve()
    if (output / 'prepared.json').exists():
        raise ValueError('Bundle already prepared; use it or choose a new output directory')
    assets_dir = output / 'assets'
    assets_dir.mkdir(parents=True, exist_ok=True)
    roles = {'bam', 'bai', 'reference', 'truth_vcf', 'truth_tbi', 'confident_bed'}
    if set(manifest['assets']) != roles or len({a['filename'] for a in manifest['assets'].values()}) != len(roles):
        raise ValueError('Asset roles or filenames are incomplete/duplicated')
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures = {k: pool.submit(download, v, assets_dir) for k, v in manifest['assets'].items()}
        assets = {k: f.result() for k, f in futures.items()}
    commands = []

    def run(args, capture=False):
        args = list(map(str, args))
        commands.append([v.replace(str(output), '<bundle>') for v in args])
        print('Running: ' + ' '.join(commands[-1][:5]), flush=True)
        return subprocess.run(args, check=True, text=True, stdout=subprocess.PIPE if capture else None).stdout

    refdir = output / 'reference'
    refdir.mkdir(exist_ok=True)
    reference = refdir / 'reference.fa'
    with gzip.open(assets['reference'], 'rb') as source, reference.open('wb') as target:
        shutil.copyfileobj(source, target, 1024 * 1024)
    run(['samtools', 'faidx', reference])
    run(['samtools', 'dict', '-o', refdir / 'reference.dict', reference])
    lengths = {f[0]: int(f[1]) for line in Path(str(reference) + '.fai').read_text().splitlines() if (f := line.split('\t'))}
    header = run(['samtools', 'view', '--no-PG', '-H', assets['bam']], True)
    contigs, samples, programs, sort_order = {}, set(), [], None
    for line in header.splitlines():
        fields = line.split('\t')
        tags = dict(f.split(':', 1) for f in fields[1:] if ':' in f)
        if fields[0] == '@SQ':
            contigs[tags['SN']] = int(tags['LN'])
        elif fields[0] == '@RG':
            samples.add(tags.get('SM'))
        elif fields[0] == '@HD':
            sort_order = tags.get('SO')
        elif fields[0] == '@PG':
            programs.append(tags)
    if samples != {manifest['sample']} or sort_order != 'coordinate':
        raise ValueError('BAM sample identity or sort order is incorrect')
    if not contigs or any(lengths.get(c) != n for c, n in contigs.items()):
        raise ValueError('BAM/reference contig dictionaries disagree')
    if run(['bcftools', 'query', '-l', assets['truth_vcf']], True).splitlines() != [manifest['sample']]:
        raise ValueError('Truth VCF sample identity is incorrect')
    region = manifest['scoring_region']
    chrom, start, end = region['chrom'], region['start'], region['end']
    if chrom not in lengths or end > lengths[chrom]:
        raise ValueError('Scoring region is outside reference bounds')
    left, right = max(0, start - manifest['padding_bases']), min(lengths[chrom], end + manifest['padding_bases'])
    interval = f'{chrom}:{left + 1}-{right}'
    (refdir / 'targets.bed').write_text(f'{chrom}\t{left}\t{right}\n')
    confident = intersect_bed(assets['confident_bed'], chrom, start, end)
    (output / 'evaluation.bed').write_text(''.join(f'{chrom}\t{s}\t{e}\n' for s, e in confident))
    run(['samtools', 'view', '--no-PG', '-b', '-o', output / 'reads.bam', assets['bam'], interval])
    run(['samtools', 'index', output / 'reads.bam'])
    run(['samtools', 'quickcheck', '-v', output / 'reads.bam'])
    run(['bcftools', 'view', '--no-version', '--regions-overlap', '1', '--regions', interval, '-Oz', '-o', output / 'truth.vcf.gz', assets['truth_vcf']])
    run(['bcftools', 'index', '--tbi', output / 'truth.vcf.gz'])
    read_count = int(run(['samtools', 'view', '-c', output / 'reads.bam'], True))
    if not read_count:
        raise ValueError('Extracted BAM contains no alignment records')
    (output / 'source_bam_header.sam').write_text(header)
    versions = {tool: run([tool, '--version'], True).strip() for tool in ('samtools', 'bcftools')}
    prepared = {'schema_version': 1, 'manifest': manifest, 'alignment_records': read_count,
        'calling_region': {'chrom': chrom, 'start': left, 'end': right}, 'evaluated_bases': sum(e - s for s, e in confident),
        'source_bam_programs': programs, 'reference_match_evidence': 'Provider assembly label and every BAM contig name/length; BAM does not contain M5 sequence hashes',
        'commands': commands, 'preparation_versions': versions, 'files': inventory(output)}
    (output / 'prepared.json').write_text(json.dumps(prepared, indent=2) + '\n')
    print(f'Prepared {manifest["sample"]}: {read_count:,} alignment records; {prepared["evaluated_bases"]:,} confident bases', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', default='benchmarks/hg002-grch38-chr20.json')
    parser.add_argument('--outdir', default='benchmark-data/hg002')
    args = parser.parse_args()
    try:
        prepare(args.manifest, args.outdir)
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(2, f'Benchmark preparation failed: {exc}\n')
