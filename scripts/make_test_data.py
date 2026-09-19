#!/usr/bin/env python3
"""Generate deterministic, entirely synthetic paired-end reads and two known SNVs."""
import argparse
import csv
import gzip
import hashlib
import json
import random
from pathlib import Path


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def write_gzip(path, text):
    with open(path, 'wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as handle:
            handle.write(text.encode('ascii'))


def generate(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20260919)
    reference = ''.join(rng.choices('ACGT', k=12000))
    (outdir / 'reference.fa').write_text('>synthetic\n' + '\n'.join(reference[i:i + 60] for i in range(0, len(reference), 60)) + '\n')
    (outdir / 'targets.bed').write_text('synthetic\t1000\t11000\n')
    truth = []
    with (outdir / 'samplesheet.csv').open('w', newline='') as sheet:
        writer = csv.writer(sheet)
        writer.writerow(['sample', 'fastq_1', 'fastq_2'])
        for sample, position, genotype in [('SYNTH_A', 4001, '0/1'), ('SYNTH_B', 8001, '1/1')]:
            ref = reference[position - 1]
            alt = next(base for base in 'ACGT' if base != ref)
            alternate = reference[:position - 1] + alt + reference[position:]
            mates = [[], []]
            for pair in range(3000):
                fragment_length = rng.randint(325, 475)
                start = rng.randrange(0, len(reference) - fragment_length + 1)
                haplotype = alternate if genotype == '1/1' or pair % 2 else reference
                sequences = [haplotype[start:start + 150], reverse_complement(haplotype[start + fragment_length - 150:start + fragment_length])]
                for mate, sequence in enumerate(sequences, 1):
                    mates[mate - 1].append(f'@{sample}_{pair}/{mate}\n{sequence}\n+\n' + 'I' * 150 + '\n')
            names = [f'{sample}_R{mate}.fastq.gz' for mate in (1, 2)]
            for name, records in zip(names, mates):
                write_gzip(outdir / name, ''.join(records))
            writer.writerow([sample, *names])
            truth.append({'sample': sample, 'chrom': 'synthetic', 'pos': position, 'ref': ref, 'alt': alt, 'gt': genotype})
    (outdir / 'expected_variants.json').write_text(json.dumps(truth, indent=2) + '\n')
    checksums = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(outdir.iterdir()) if p.name != 'provenance.json' and p.is_file()}
    (outdir / 'provenance.json').write_text(json.dumps({'kind': 'synthetic', 'seed': 20260919, 'read_pairs_per_sample': 3000,
        'purpose': 'Software integration testing only; not human accuracy benchmarking', 'sha256': checksums}, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', default='tests/data')
    generate(parser.parse_args().outdir)
