#!/usr/bin/env python3
"""Stream both mates: validate structure, matching IDs, and equal read counts."""
import argparse
import gzip
import json
import re
from contextlib import ExitStack
from itertools import zip_longest
from pathlib import Path


def open_fastq(path):
    path = Path(path)
    with path.open('rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    if compressed != (path.suffix == '.gz'):
        raise ValueError(f'{path}: gzip content must have a .gz extension, and vice versa')
    return (gzip.open if compressed else open)(path, 'rt', encoding='ascii')


def records(handle, label, mate):
    number = 0
    while True:
        header = handle.readline()
        if not header:
            return
        number += 1
        sequence = handle.readline().rstrip('\r\n')
        plus = handle.readline().rstrip('\r\n')
        quality = handle.readline().rstrip('\r\n')
        fields = header.strip().split()
        if not fields or not fields[0].startswith('@') or len(fields[0]) < 2:
            raise ValueError(f'{label}: record {number}: invalid FASTQ header')
        if not sequence or not re.fullmatch('[ACGTNacgtn]+', sequence):
            raise ValueError(f'{label}: record {number}: empty or invalid DNA sequence')
        if not plus.startswith('+') or len(sequence) != len(quality):
            raise ValueError(f'{label}: record {number}: truncated record or sequence/quality length mismatch')
        if any(ord(c) < 33 or ord(c) > 126 for c in quality):
            raise ValueError(f'{label}: record {number}: invalid quality encoding')
        read_id = fields[0][1:]
        if read_id.endswith(('/1', '/2')):
            if read_id[-1] != str(mate):
                raise ValueError(f'{label}: record {number}: wrong mate orientation')
            read_id = read_id[:-2]
        if len(fields) > 1 and fields[1].startswith(('1:', '2:')) and fields[1][0] != str(mate):
            raise ValueError(f'{label}: record {number}: wrong CASAVA mate orientation')
        yield read_id


def validate_pair(read1, read2):
    count = 0
    with ExitStack() as stack:
        r1 = stack.enter_context(open_fastq(read1))
        r2 = stack.enter_context(open_fastq(read2))
        for count, pair in enumerate(zip_longest(records(r1, 'R1', 1), records(r2, 'R2', 2)), 1):
            if None in pair:
                raise ValueError(f'Unequal FASTQ record counts at pair {count}')
            if pair[0] != pair[1]:
                raise ValueError(f'Mismatched FASTQ IDs at pair {count}: {pair[0]} vs {pair[1]}')
    if not count:
        raise ValueError('FASTQ files contain no reads')
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--read1', required=True)
    parser.add_argument('--read2', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    try:
        count = validate_pair(args.read1, args.read2)
    except (ValueError, OSError, EOFError, UnicodeError) as exc:
        parser.exit(2, f'FASTQ validation failed for {args.sample}: {exc}\n')
    Path(args.output).write_text(json.dumps({'sample': args.sample, 'read_pairs': count, 'validated': True}, indent=2) + '\n')


if __name__ == '__main__':
    main()
