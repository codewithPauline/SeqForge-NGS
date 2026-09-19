#!/usr/bin/env python3
"""Exercise Nextflow's actual samplesheet validation before expensive tasks launch."""
import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def check(nextflow):
    root = Path(__file__).resolve().parents[1]
    fixture = root / 'tests/data'
    row = ['DUPLICATE', str(fixture / 'SYNTH_A_R1.fastq.gz'), str(fixture / 'SYNTH_A_R2.fastq.gz')]
    cases = [
        ('duplicate', [row, row], 'Duplicate sample ID'),
        ('reuse', [row, ['OTHER', *row[1:]]], 'FASTQ reused'),
        ('empty', [], 'Samplesheet contains no samples'),
        ('unsafe', [['bad;id', *row[1:]]], 'Invalid sample ID'),
    ]
    for name, rows, message in cases:
        with tempfile.TemporaryDirectory(prefix='seqforge-validation-') as temp:
            temp = Path(temp)
            sheet = temp / 'samples.csv'
            with sheet.open('w', newline='') as handle:
                writer = csv.writer(handle)
                writer.writerow(['sample', 'fastq_1', 'fastq_2'])
                writer.writerows(rows)
            result = subprocess.run([nextflow, 'run', str(root / 'main.nf'),
                '--input', str(sheet), '--fasta', str(fixture / 'reference.fa'),
                '--outdir', str(temp / 'results'), '-work-dir', str(temp / 'work'), '-ansi-log', 'false'],
                cwd=temp, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
            assert result.returncode != 0 and message in result.stdout, f'{name}: {result.stdout}'
            assert 'Submitted process' not in result.stdout, f'{name}: invalid input launched a task'
            print(f'Samplesheet rejection verified: {name}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nextflow', default='nextflow')
    check(parser.parse_args().nextflow)
