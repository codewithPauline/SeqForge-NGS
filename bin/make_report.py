#!/usr/bin/env python3
"""Summarize measured QC. Variant counts are not accuracy estimates."""
import argparse
import html
import json
from pathlib import Path


def samtools_stats(path):
    values = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith('SN\t'):
            fields = line.split('\t')
            values[fields[1].rstrip(':')] = float(fields[2])
    required = {'raw total sequences', 'reads mapped', 'reads duplicated'}
    if not required <= values.keys():
        raise ValueError(f'Missing required alignment metrics in {path}')
    return values


def variant_stats(path):
    values = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith('SN\t'):
            fields = line.split('\t')
            values[fields[2].rstrip(':')] = int(fields[3])
    required = {'number of records', 'number of SNPs', 'number of indels'}
    if not required <= values.keys():
        raise ValueError(f'Missing required variant metrics in {path}')
    return values


def summarize(metrics):
    samples = []
    for file in sorted(Path(metrics).glob('*.fastp.json')):
        sample = file.name.removesuffix('.fastp.json')
        fastp = json.loads(file.read_text())['summary']
        alignment = samtools_stats(Path(metrics) / f'{sample}.samtools.stats')
        variants = variant_stats(Path(metrics) / f'{sample}.bcftools.stats')
        total = alignment['raw total sequences']
        samples.append({'sample': sample,
            'raw_reads': fastp['before_filtering']['total_reads'],
            'retained_reads': fastp['after_filtering']['total_reads'],
            'q30_rate_after_filtering': fastp['after_filtering']['q30_rate'],
            'mapped_percent': round(100 * alignment['reads mapped'] / total, 3) if total else None,
            'duplicate_percent': round(100 * alignment['reads duplicated'] / total, 3) if total else None,
            'pass_records': variants['number of records'],
            'pass_snps': variants['number of SNPs'], 'pass_indels': variants['number of indels']})
    if not samples:
        raise ValueError('Cannot generate a report without sample metrics')
    return samples


def render_report(report):
    rows = []
    for sample in report['samples']:
        mapped = sample['mapped_percent']
        mapped_label = f'{mapped:.2f}%' if mapped is not None else 'N/A'
        rows.append('<tr>' + ''.join(f'<td>{html.escape(str(value))}</td>' for value in (
            sample['sample'], f"{sample['raw_reads']:,}", f"{sample['retained_reads']:,}",
            mapped_label, sample['pass_snps'], sample['pass_indels'])) + '</tr>')
    metadata = html.escape(json.dumps({k: v for k, v in report.items() if k != 'samples'}, indent=2))
    return '''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SeqForge-NGS | Run summary</title>
<style>
body{font-family:system-ui,sans-serif;margin:0;background:#f5f7fa;color:#172c3c}
main{max-width:1100px;margin:48px auto;padding:0 28px}h1{font-size:42px;margin-bottom:8px}
.eyebrow{color:#087e83;font-weight:700;letter-spacing:.1em;font-size:13px}
.card{background:white;padding:24px;border:1px solid #dce5eb;border-radius:12px;margin:24px 0;overflow:auto}
table{width:100%;border-collapse:collapse;text-align:left}th,td{padding:14px 12px;border-bottom:1px solid #e4eaf0}
th{font-size:12px;text-transform:uppercase;color:#526b80}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}
p{line-height:1.6;color:#526b80}a{color:#087e83}footer{font-size:13px;margin:32px 0}
</style><main><div class="eyebrow">GERMLINE VARIANT DISCOVERY</div><h1>SeqForge-NGS</h1>
<p>Measured read quality, alignment, and PASS variant counts from this run.</p>
<div class="card"><table><thead><tr><th>Sample</th><th>Raw reads</th><th>Retained reads</th>
<th>Mapped</th><th>PASS SNPs</th><th>PASS indels</th></tr></thead><tbody>''' + ''.join(rows) + '''</tbody></table></div>
<p>Counts describe pipeline output; they do not establish precision, recall, or biological validity.
The baseline filter uses QUAL and per-sample DP. Review <a href="multiqc_report.html">MultiQC</a>
and the raw/filtered VCFs before interpreting results.</p>
<div class="card"><h2>Run provenance</h2><pre>''' + metadata + '''</pre></div>
<footer>SeqForge-NGS · Pauline Owusu-Ansah · Research software in development</footer></main></html>'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('metrics', 'versions', 'reference', 'pipeline-version', 'nextflow-version'):
        parser.add_argument(f'--{name}', required=True)
    parser.add_argument('--min-qual', type=float, required=True)
    parser.add_argument('--min-depth', type=int, required=True)
    args = parser.parse_args()
    report = {'pipeline_version': args.pipeline_version, 'nextflow_version': args.nextflow_version,
        'reference': json.loads(Path(args.reference).read_text()),
        'caller': 'GATK HaplotypeCaller; diploid single-sample VCF',
        'filter': {'name': 'LowQualDepth', 'min_qual': args.min_qual, 'min_sample_depth': args.min_depth,
                   'missing_qual_or_depth': 'filtered'},
        'samples': summarize(args.metrics)}
    Path('seqforge_summary.json').write_text(json.dumps(report, indent=2) + '\n')
    Path('seqforge_summary.html').write_text(render_report(report))
    versions = '\n'.join(f'[{p.name}]\n{p.read_text().strip()}\n' for p in sorted(Path(args.versions).glob('*.txt')))
    if not versions:
        raise ValueError('No tool version records found')
    Path('software_versions.txt').write_text(versions)


if __name__ == '__main__':
    main()
