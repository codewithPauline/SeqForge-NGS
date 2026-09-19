#!/usr/bin/env python3
"""Check real workflow outputs against the synthetic truth, including genotypes."""
import argparse
import csv
import gzip
import json
import subprocess
from pathlib import Path


def check(results, fixture, resumed=False):
    results, fixture = Path(results), Path(fixture)
    truth = json.loads((fixture / 'expected_variants.json').read_text())
    summary = json.loads((results / 'reports/seqforge_summary.json').read_text())
    assert {s['sample'] for s in summary['samples']} == {v['sample'] for v in truth}
    for expected in truth:
        sample = expected['sample']
        vcf = results / f'variants/filtered/{sample}.pass.vcf.gz'
        records = []
        with gzip.open(vcf, 'rt') as handle:
            for line in handle:
                if line.startswith('#CHROM'):
                    assert line.rstrip().split('\t')[9:] == [sample], 'VCF sample identity was lost'
                elif not line.startswith('#'):
                    fields = line.rstrip().split('\t')
                    values = dict(zip(fields[8].split(':'), fields[9].split(':')))
                    records.append((fields[0], int(fields[1]), fields[3], fields[4],
                        '/'.join(sorted(values['GT'].replace('|', '/').split('/')))))
        key = tuple(expected[k] for k in ('chrom', 'pos', 'ref', 'alt', 'gt'))
        assert records == [key], f'{sample}: expected {key}, observed {records}'
        bam = results / f'alignment/{sample}.bam'
        subprocess.run(['samtools', 'quickcheck', '-v', str(bam)], check=True)
        header = subprocess.check_output(['samtools', 'view', '-H', str(bam)], text=True)
        assert f'SM:{sample}' in header and 'SO:coordinate' in header
        assert int(subprocess.check_output(['samtools', 'view', '-c', str(bam), 'synthetic:3900-8100'])) > 0
        assert (results / f'alignment/{sample}.bam.bai').is_file()
        assert Path(str(vcf) + '.tbi').is_file()
        qc = next(s for s in summary['samples'] if s['sample'] == sample)
        assert qc['raw_reads'] == 6000
        assert qc['mapped_percent'] > 95
        assert qc['pass_records'] == 1 and qc['pass_snps'] == 1 and qc['pass_indels'] == 0
    for name in ('reports/multiqc_report.html', 'reports/seqforge_summary.html',
                 'pipeline_info/software_versions.txt', 'pipeline_info/reference.json',
                 'pipeline_info/execution_report.html', 'pipeline_info/timeline.html'):
        assert (results / name).stat().st_size > 0, f'Missing output: {name}'
    if resumed:
        with (results / 'pipeline_info/trace.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        assert len(rows) == 17, f'Expected 17 executed tasks, observed {len(rows)}'
        assert {row['status'] for row in rows} == {'CACHED'}, 'Resume reran unchanged tasks'
    print(f'Integration checks passed: {len(truth)} samples, expected variants/genotypes, BAM identities, QC, reports' + ('; all tasks cached' if resumed else ''))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', default='results-test')
    parser.add_argument('--fixture', default='tests/data')
    parser.add_argument('--resumed', action='store_true')
    args = parser.parse_args()
    check(args.results, args.fixture, args.resumed)
