#!/usr/bin/env python3
"""Validate FASTA and BED compatibility without retaining a genome in memory."""
import argparse
import hashlib
import json
import re
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def fasta_lengths(path):
    contigs = {}
    current = None
    with open(path, encoding='ascii') as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                fields = line[1:].split()
                if not fields or not re.fullmatch(r'[A-Za-z0-9_.:-]+', fields[0]):
                    raise ValueError(f'Invalid FASTA contig name at line {number}')
                current = fields[0]
                if current in contigs:
                    raise ValueError(f'Duplicate FASTA contig: {current}')
                contigs[current] = 0
            else:
                if current is None or not re.fullmatch('[ACGTNRYKMSWBDHVacgtnrykmswbdhv]+', line):
                    raise ValueError(f'Invalid FASTA sequence at line {number}')
                contigs[current] += len(line)
    if not contigs or any(length == 0 for length in contigs.values()):
        raise ValueError('Reference contains no sequence or an empty contig')
    return contigs


def bed_intervals(path, contigs):
    intervals = []
    order = {name: i for i, name in enumerate(contigs)}
    previous = None
    with open(path) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith('#'):
                continue
            fields = raw.split()
            if len(fields) < 3:
                raise ValueError(f'BED line {number}: expected at least 3 columns')
            chrom = fields[0]
            try:
                start, end = map(int, fields[1:3])
            except ValueError as exc:
                raise ValueError(f'BED line {number}: coordinates must be integers') from exc
            if chrom not in contigs:
                raise ValueError(f'BED contig {chrom} not found in reference (check chr naming and assembly)')
            if not 0 <= start < end <= contigs[chrom]:
                raise ValueError(f'BED line {number}: interval is empty or outside reference bounds')
            key = (order[chrom], start, end)
            if previous and (key < previous or (key[0] == previous[0] and start < previous[2])):
                raise ValueError('BED intervals must be non-overlapping and sorted in FASTA contig order')
            previous = key
            intervals.append((chrom, start, end))
    if not intervals:
        raise ValueError('BED contains no intervals')
    return intervals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--intervals')
    parser.add_argument('--bed-out', required=True)
    parser.add_argument('--metadata', required=True)
    args = parser.parse_args()
    try:
        contigs = fasta_lengths(args.fasta)
        intervals = bed_intervals(args.intervals, contigs) if args.intervals else [(c, 0, n) for c, n in contigs.items()]
    except (ValueError, OSError, UnicodeError) as exc:
        parser.exit(2, f'Reference validation failed: {exc}\n')
    Path(args.bed_out).write_text(''.join(f'{c}\t{s}\t{e}\n' for c, s, e in intervals))
    metadata = {'reference_sha256': sha256(args.fasta), 'contigs': contigs,
                'scope': 'provided_intervals' if args.intervals else 'whole_reference',
                'target_bases': sum(e - s for _, s, e in intervals),
                'intervals_sha256': sha256(args.bed_out)}
    Path(args.metadata).write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    main()
