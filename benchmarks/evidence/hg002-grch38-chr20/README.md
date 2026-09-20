# HG002 regional calling benchmark

Measured on 2026-09-20 using the committed resource manifest and the source-file hashes in `source_snapshot.json`.

**Scope:** calling and baseline filtering from public HG002 BWA/Picard alignments. SeqForge's FASTQ QC, trimming, BWA-MEM2 alignment, and duplicate marking are not evaluated here.

- Reference: full GRCh38 no-alt analysis set, GCA_000001405.15.
- Scoring interval: chr20:10,000,001–11,000,000 (1-based inclusive).
- Confident territory: 977,839 bases after intersection with GIAB v4.2.1.
- Calling interval: chr20:9,990,001–11,010,000, providing 10 kb flanks.
- Input: 267,154 existing alignment records in the padded interval; source labeled NovaSeq PCR-free 35x.
- Caller: GATK 4.6.1.0 HaplotypeCaller, diploid, no BQSR.
- Filter: preset QUAL ≥ 30 and sample DP ≥ 10; missing values filtered.
- Evaluator: RTG vcfeval 3.12.1 with genotype matching and transborder matches.

| Type | Precision | Recall | F1 | TP baseline | TP call | FP | FN |
|---|---|---|---|---|---|---|---|
| SNP | 99.632% | 99.852% | 99.742% | 1352 | 1353 | 5 | 2 |
| INDEL | 96.356% | 96.581% | 96.469% | 226 | 238 | 9 | 8 |

TP baseline and TP call use RTG's distinct truth/call representations. Precision uses TP call; recall uses TP baseline. Variant-type baseline counts are RTG-rescaled counts, not naïve VCF line counts. Metrics use the final all-PASS ROC point; the automatically printed best-F1 row in `rtg_summary.txt` is not used for reporting.

This is one prespecified region in one sample, not a genome-wide estimate or independent clinical validation. No filter thresholds were tuned to improve these results. Indel errors remain visible and should be investigated on separate development regions before evaluating changes on held-out regions.

## Evidence

- `benchmark_summary.json` and `.html`: measured metrics, limitations, versions, and source hashes.
- `prepared.json`: pinned source URLs, provider MD5 checksums, preparation commands, and derived-file SHA-256 values.
- `weighted_roc.tsv`, `snp_roc.tsv`, `indel_roc.tsv`: original RTG output, decompressed without editing.
- `fp.vcf`, `fn.vcf`: original RTG false-positive and false-negative records, decompressed without editing.
- `rtg_summary.txt`: original summary including the all-calls (`None`) row.

Large FASTA/BAM resources are reproducibly downloaded, not committed. Follow [the benchmark instructions](../../../docs/benchmark.md) to regenerate the run. The provider's chromosome 20 BAM already reflects selection by alignment position, which limits extrapolation to other regions and raw-read workflows.
