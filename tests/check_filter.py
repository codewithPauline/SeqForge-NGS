#!/usr/bin/env python3
"""Verify inclusive thresholds and missing values using the actual filter script."""
import gzip
import subprocess
import tempfile
from pathlib import Path


def check():
    script = Path(__file__).resolve().parents[1] / 'bin/filter_variants.sh'
    with tempfile.TemporaryDirectory(prefix='seqforge-filter-') as temp:
        root = Path(temp)
        reference = root / 'reference.fa'
        reference.write_text('>chr1\n' + 'A' * 100 + '\n')
        subprocess.run(['samtools', 'faidx', str(reference)], check=True)
        header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTEST\n')
        cases = [('30', '10'), ('29', '20'), ('60', '9'), ('.', '20'), ('60', '.'), ('60', '20')]
        vcf = root / 'input.vcf'
        vcf.write_text(header + ''.join(f'chr1\t{i}\t.\tA\tC\t{qual}\t.\t.\tGT:DP\t0/1:{depth}\n'
            for i, (qual, depth) in enumerate(cases, 1)))
        subprocess.run(['bash', str(script), str(reference), str(vcf), 'TEST', '30', '10'], cwd=root, check=True)
        with gzip.open(root / 'TEST.filtered.vcf.gz', 'rt') as handle:
            filters = {int(line.split('\t')[1]): line.split('\t')[6] for line in handle if not line.startswith('#')}
        assert filters == {1: 'PASS', 2: 'LowQualDepth', 3: 'LowQualDepth', 4: 'LowQualDepth', 5: 'LowQualDepth', 6: 'PASS'}, filters
        with gzip.open(root / 'TEST.pass.vcf.gz', 'rt') as handle:
            positions = [int(line.split('\t')[1]) for line in handle if not line.startswith('#')]
        assert positions == [1, 6], positions
    print('Filter checks passed: boundary QUAL/DP, low values, missing values, PASS extraction')


if __name__ == '__main__':
    check()
