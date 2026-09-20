#!/usr/bin/env python3
"""Report the all-PASS operating point from RTG, without tuning thresholds on truth."""
import argparse
import gzip
import html
import json
import math
from pathlib import Path
from benchmark_io import checksum


def read_roc(path):
    baseline = calls = None
    rows, header = [], None
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#total baseline variants:'):
                baseline = float(line.split(':', 1)[1])
            elif line.startswith('#total call variants:'):
                calls = float(line.split(':', 1)[1])
            elif line.startswith('#score\t'):
                header = line.strip().lstrip('#').split('\t')
            elif line.strip() and not line.startswith('#'):
                fields = line.strip().split('\t')
                if not header or len(fields) != len(header):
                    raise ValueError(f'Invalid ROC header/row: {path}')
                rows.append(dict(zip(header, fields)))
    if baseline is None or calls is None or header is None:
        raise ValueError(f'Incomplete ROC file: {path}')
    if rows:
        # Descending QUAL: the final operating point includes all scored calls.
        row = rows[-1]
        tp_baseline, tp_call, fp = (float(row[k]) for k in ('true_positives_baseline', 'true_positives_call', 'false_positives'))
        fn = float(row.get('false_negatives', baseline - tp_baseline))
        if not math.isclose(tp_call + fp, calls, abs_tol=0.011) or not math.isclose(tp_baseline + fn, baseline, abs_tol=0.011):
            raise ValueError(f'Final ROC point does not account for all calls/truth: {path}')
        threshold = row['score']
    else:
        if calls != 0:
            raise ValueError(f'ROC contains calls but no operating points: {path}')
        tp_baseline = tp_call = fp = 0.0
        fn, threshold = baseline, None
    precision = tp_call / (tp_call + fp) if tp_call + fp else None
    recall = tp_baseline / (tp_baseline + fn) if tp_baseline + fn else None
    f1 = None if precision is None or recall is None else (2 * precision * recall / (precision + recall) if precision + recall else 0.0)
    return {'truth_total': baseline, 'query_total': calls, 'tp_baseline': tp_baseline, 'tp_call': tp_call, 'fp': fp, 'fn': fn,
            'precision': precision, 'recall': recall, 'f1': f1, 'lowest_included_qual': threshold, 'operating_point': 'all_scored_PASS_calls'}


def build_report(comparison, prepared, query, versions, source):
    comparison = Path(comparison)
    provenance = json.loads(Path(prepared).read_text())
    manifest = provenance['manifest']
    return {'schema_version': 1, 'benchmark_id': manifest['id'], 'sample': manifest['sample'],
        'scope': manifest['analysis_scope'], 'assembly': manifest['assembly'], 'truth_release': manifest['truth_release'],
        'scoring_region': manifest['scoring_region'], 'calling_region': provenance['calling_region'],
        'evaluated_bases': provenance['evaluated_bases'], 'alignment_records': provenance['alignment_records'],
        'filter': {'min_qual': manifest['min_qual'], 'min_sample_depth': manifest['min_depth']},
        'engine': 'RTG vcfeval 3.12.1; diploid genotype matching; transborder matches allowed',
        'metric_definition': 'RTG weighted ROC counts, with baseline rescaling for variant-type subsets. Final all-PASS point; no best-F1 threshold selection.',
        'metrics': {label: read_roc(comparison / name) for label, name in [('ALL', 'weighted_roc.tsv.gz'), ('SNP', 'snp_roc.tsv.gz'), ('INDEL', 'indel_roc.tsv.gz')]},
        'provenance': {'prepared_json_sha256': checksum(prepared), 'query_vcf_sha256': checksum(query),
            'source_files': json.loads(Path(source).read_text()),
            'software_versions': {p.name: p.read_text().strip() for p in sorted(Path(versions).glob('*.txt'))}},
        'limitations': ['One selected genomic region and one sample; not a genome-wide accuracy estimate.',
            'Starts from existing alignments; SeqForge FASTQ QC, trimming, alignment, and duplicate marking are not evaluated.',
            'Baseline QUAL/depth filter; no BQSR or threshold tuning in this case study.',
            'Variant-type counts follow RTG classification and baseline rescaling and may not partition complex variants.']}


def render(report):
    def metric(value):
        return 'N/A' if value is None else f'{100 * value:.3f}%'
    rows = ''.join('<tr>' + ''.join(f'<td>{html.escape(str(x))}</td>' for x in (label, v['tp_baseline'], v['tp_call'], v['fp'], v['fn'], metric(v['precision']), metric(v['recall']), metric(v['f1']))) + '</tr>' for label, v in report['metrics'].items())
    region = report['scoring_region']
    details = html.escape(f"{report['sample']} · {report['assembly']} · {region['chrom']}:{region['start'] + 1:,}–{region['end']:,}")
    limitations = ''.join(f'<li>{html.escape(v)}</li>' for v in report['limitations'])
    return '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>SeqForge-NGS benchmark</title><style>body{font-family:system-ui,sans-serif;margin:0;background:#f5f7fa;color:#193245}
main{max-width:1180px;margin:48px auto;padding:0 28px}h1{font-size:40px}p,li{line-height:1.65}.card{background:white;padding:24px;border:1px solid #d7e1e8;border-radius:12px;overflow:auto}
table{width:100%;border-collapse:collapse;text-align:left}th,td{padding:12px;border-bottom:1px solid #e1e8ee}th{font-size:12px;color:#526c80}.label{color:#087e83;font-weight:700}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}</style>
<main><p class="label">REGIONAL CALLING AND FILTERING BENCHMARK</p><h1>SeqForge-NGS</h1><p>''' + details + '''</p>
<div class="card"><table><thead><tr><th>Type</th><th>TP baseline</th><th>TP call</th><th>FP</th><th>FN</th><th>Precision</th><th>Recall</th><th>F1</th></tr></thead><tbody>''' + rows + '''</tbody></table></div>
<p>All scored PASS calls at the preset filter thresholds. N/A means the metric has no denominator.</p><h2>Scope and limitations</h2><ul>''' + limitations + '''</ul><details><summary>Provenance and complete metrics</summary><pre>''' + html.escape(json.dumps(report, indent=2)) + '</pre></details></main></html>'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ('comparison', 'prepared', 'query', 'versions', 'source'):
        parser.add_argument('--' + argument, required=True)
    args = parser.parse_args()
    report = build_report(args.comparison, args.prepared, args.query, args.versions, args.source)
    Path('benchmark_summary.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    Path('benchmark_summary.html').write_text(render(report))
