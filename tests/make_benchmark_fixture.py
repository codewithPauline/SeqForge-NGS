#!/usr/bin/env python3
"""Prepare a tiny caller benchmark with a known SNV and an intentionally missed indel."""
import json
import shutil
import subprocess
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from benchmark_io import inventory


def prepare():
    bundle = Path('benchmark-data/ci')
    refdir = bundle / 'reference'
    refdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile('tests/data/reference.fa', refdir / 'reference.fa')
    shutil.copyfile('tests/data/targets.bed', refdir / 'targets.bed')
    shutil.copyfile('tests/data/targets.bed', bundle / 'evaluation.bed')
    shutil.copyfile('results-test/alignment/SYNTH_A.bam', bundle / 'reads.bam')
    shutil.copyfile('results-test/alignment/SYNTH_A.bam.bai', bundle / 'reads.bam.bai')
    subprocess.run(['samtools', 'faidx', str(refdir / 'reference.fa')], check=True)
    subprocess.run(['samtools', 'dict', '-o', str(refdir / 'reference.dict'), str(refdir / 'reference.fa')], check=True)
    reference = ''.join(line for line in (refdir / 'reference.fa').read_text().splitlines() if not line.startswith('>'))
    snv = json.loads(Path('tests/data/expected_variants.json').read_text())[0]
    records = [(snv['pos'], snv['ref'], snv['alt']), (7001, reference[7000:7002], reference[7000])]
    truth = bundle / 'truth.vcf'
    truth.write_text('##fileformat=VCFv4.2\n##contig=<ID=synthetic,length=12000>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH_A\n' +
        ''.join(f'synthetic\t{p}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT\t0/1\n' for p, ref, alt in records))
    subprocess.run(['bcftools', 'view', '-Oz', '-o', str(bundle / 'truth.vcf.gz'), str(truth)], check=True)
    subprocess.run(['bcftools', 'index', '--tbi', str(bundle / 'truth.vcf.gz')], check=True)
    manifest = {'schema_version': 1, 'id': 'synthetic-ci', 'sample': 'SYNTH_A',
        'analysis_scope': 'synthetic_software_test_only', 'assembly': 'synthetic', 'truth_release': 'generated-ci-fixture',
        'scoring_region': {'chrom': 'synthetic', 'start': 1000, 'end': 11000, 'coordinate_system': '0-based half-open'},
        'padding_bases': 0, 'min_qual': 30, 'min_depth': 10}
    prepared = {'schema_version': 1, 'manifest': manifest, 'calling_region': manifest['scoring_region'],
                'alignment_records': 6000, 'evaluated_bases': 10000, 'files': inventory(bundle)}
    (bundle / 'prepared.json').write_text(json.dumps(prepared, indent=2) + '\n')


if __name__ == '__main__':
    prepare()
