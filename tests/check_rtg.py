#!/usr/bin/env python3
"""Real evaluator regression: indel representation, genotype mismatch, and BED masking."""
import os
import subprocess
import sys
import tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from benchmark_report import read_roc


def check():
    with tempfile.TemporaryDirectory(prefix='seqforge-rtg-') as tmp:
        root = Path(tmp)
        sequence = 'ACGT' * 20 + 'A' * 20 + 'CGT' * 70
        reference = root / 'reference.fa'
        reference.write_text('>chr1\n' + sequence + '\n')
        bed = root / 'evaluation.bed'
        bed.write_text('chr1\t0\t150\n')
        header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=310>\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTEST\n')
        truth = [(30, sequence[29], 'G', '0/1'), (85, 'AA', 'A', '0/1')]
        query = [(30, sequence[29], 'G', '1/1'), (90, 'AA', 'A', '0/1'), (200, sequence[199], 'A', '0/1')]
        for name, records in [('truth', truth), ('query', query)]:
            path = root / f'{name}.vcf'
            path.write_text(header + ''.join(f'chr1\t{p}\t.\t{r}\t{a}\t60\tPASS\t.\tGT\t{gt}\n' for p, r, a, gt in records))
            subprocess.run(['bcftools', 'view', '-Oz', '-o', str(path) + '.gz', str(path)], check=True)
            subprocess.run(['bcftools', 'index', '--tbi', str(path) + '.gz'], check=True)
        environment = dict(os.environ, RTG_MEM='1g')
        subprocess.run(['rtg', 'format', '-o', str(root / 'ref.sdf'), str(reference)], env=environment, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(['rtg', 'vcfeval', '-b', str(root / 'truth.vcf.gz'), '-c', str(root / 'query.vcf.gz'), '-t', str(root / 'ref.sdf'), '-e', str(bed), '-o', str(root / 'comparison'), '--sample', 'TEST', '--vcf-score-field', 'QUAL', '--roc-subset', 'snp,indel', '-T', '2'], env=environment, check=True)
        results = read_roc(root / 'comparison/weighted_roc.tsv.gz')
        assert (results['tp_baseline'], results['tp_call'], results['fp'], results['fn']) == (1, 1, 1, 1), results
        assert read_roc(root / 'comparison/indel_roc.tsv.gz')['recall'] == 1
        assert read_roc(root / 'comparison/snp_roc.tsv.gz')['recall'] == 0
    print('RTG checks passed: equivalent indels match, wrong genotypes fail, outside-BED calls excluded')


if __name__ == '__main__':
    check()
