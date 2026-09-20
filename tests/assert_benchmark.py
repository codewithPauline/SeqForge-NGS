#!/usr/bin/env python3
"""Check the benchmark workflow, including zero-call indel metrics and resume."""
import argparse
import csv
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--resumed', action='store_true')
args = parser.parse_args()
results = Path('results-benchmark-test')
report = json.loads((results / 'benchmark/benchmark_summary.json').read_text())
assert report['sample'] == 'SYNTH_A' and report['evaluated_bases'] == 10000
assert report['scope'] == 'synthetic_software_test_only'
snps, indels = report['metrics']['SNP'], report['metrics']['INDEL']
assert (snps['tp_baseline'], snps['tp_call'], snps['fp'], snps['fn']) == (1, 1, 0, 0), snps
assert (indels['tp_baseline'], indels['tp_call'], indels['fp'], indels['fn']) == (0, 0, 0, 1), indels
assert indels['precision'] is None and indels['recall'] == 0
for name in ('benchmark_summary.html', 'prepared.json', 'source_snapshot.json', 'comparison/summary.txt'):
    assert (results / 'benchmark' / name).stat().st_size > 0
if args.resumed:
    with (results / 'pipeline_info/trace.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 6 and {r['status'] for r in rows} == {'CACHED'}, rows
print('Benchmark integration passed: expected SNP, missed indel, honest empty-call metrics, provenance' + ('; six cached tasks' if args.resumed else ''))
