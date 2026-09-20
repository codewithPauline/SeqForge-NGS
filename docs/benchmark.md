# Reproduce the HG002 benchmark

This benchmark evaluates **GATK calling and baseline filtering from existing public alignments**. The ordinary `main.nf` FASTQ workflow is separate. Its input validation, trimming, alignment, and duplicate marking are tested by the synthetic integration suite, but are not assessed by the human accuracy scores below.

See [measured results and evidence](../benchmarks/evidence/hg002-grch38-chr20/README.md).

## Fixed resources and scope

The [resource manifest](../benchmarks/hg002-grch38-chr20.json) pins public Google Cloud Storage object generations, byte sizes, and provider MD5 checksums for:

- HG002 NovaSeq PCR-free, nominal 35x, chromosome 20 BAM and index. Its header records BWA 0.7.17 and Picard/GATK 4.1.2.0 duplicate marking.
- Full GRCh38 no-alt analysis reference, GCA_000001405.15.
- HG002 GIAB v4.2.1 truth VCF, index, and matching confident-region BED.

The [DeepVariant public case-study documentation](https://github.com/google/deepvariant/blob/r1.8/docs/deepvariant-case-study.md) describes this resource family; its worked example uses HG003. SeqForge explicitly selects HG002 objects and verifies that both BAM and truth VCF identify HG002. [NIST's HG002 release](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/) is the upstream benchmark source. This pipeline runs GATK, not DeepVariant.

Scoring covers chr20:10,000,001–11,000,000 (1-based inclusive), intersected with confident regions. A 10 kb flank is retained for calling and haplotype comparison. Reads are extracted by alignment overlap; original alignments, read groups, and duplicate flags are preserved. No reads are realigned to a cropped reference.

## Run with Docker

Use Linux x86-64, Nextflow 24.10.6, Java 17, approximately 16 GB RAM, and at least 20 GB free disk including the image and working files. The pinned downloads total about 2.15 GB; extraction includes a full 3.14 GB FASTA. These are approximate preparation requirements, not WGS performance claims.

```bash
docker build -t seqforge-ngs:0.1.0-dev .
docker run --rm -v "$PWD:$PWD" -w "$PWD" \
  --entrypoint /bin/bash seqforge-ngs:0.1.0-dev \
  -c 'python scripts/prepare_benchmark.py'

nextflow run benchmark.nf -profile docker \
  --benchmark_data benchmark-data/hg002 \
  --outdir results-hg002 --max_cpus 4 --calling_memory '6 GB'
```

With the locked environment installed locally, run `python scripts/prepare_benchmark.py` directly and omit `-profile docker`. Thresholds come from the prepared manifest, not `--min_qual` or `--min_depth` CLI options. Change the manifest and prepare a separate bundle for an explicitly different experiment; retain the original benchmark as a fixed baseline.

The preparation script checks every downloaded object before using it, verifies BAM sorting and sample identity, compares all BAM contig names/lengths with the reference, checks truth sample identity, and records checksums for the derived files. The source BAM lacks M5 sequence hashes, so reference matching relies on provider metadata and dictionary agreement rather than claiming sequence-hash proof from the BAM.

Already downloaded, valid assets can be reused after an interrupted preparation. A completed `prepared.json` is not overwritten. A corrupt cached download fails clearly and must be removed before retrying. The workflow verifies derived SHA-256 values before running expensive tasks.

## Interpretation

The benchmark reuses the same `HAPLOTYPECALLER` and `FILTER_VARIANTS` modules as `main.nf`. RTG vcfeval compares diploid genotypes while accounting for equivalent haplotype representations. `--evaluation-regions` permits transborder matches; reference and truth context remain available outside scoring boundaries.

Reports use the **last descending-QUAL ROC point**, including every scored PASS call at the preset QUAL/depth filters. RTG's default console summary also prints a best-F1 threshold; that optimized row is not used for reported metrics. The parser checks that the selected ROC point accounts for all calls and baseline variants.

Precision = TP-call / (TP-call + FP). Recall = TP-baseline / (TP-baseline + FN). RTG can match different numbers of truth and call records for the same haplotype, so these TP counts remain distinct. SNP/INDEL subsets use RTG's baseline rescaling. Undefined denominators remain null/N/A rather than becoming 100%.

This single-region result does not establish genome-wide accuracy, performance on other ancestries or platforms, or a complete FASTQ-to-VCF validation. It contains no BQSR, VQSR, clinical interpretation, or functional annotation. Evaluate future tuning on separate development regions and reserve held-out regions for evaluation.

## Outputs and automated checks

`results-hg002/benchmark/` contains HTML/JSON summaries, prepared-data provenance, source-file hashes, and the original RTG comparison directory. Raw/filtered/PASS VCFs and Nextflow execution reports use the existing output layout.

CI runs a small benchmark using a known synthetic SNP and an intentionally absent indel, checking correct true-positive/false-negative accounting and honest empty-call metrics. An independent evaluator test checks that shifted but equivalent deletions match, wrong genotypes fail, and calls outside confident regions are excluded. CI also verifies that all six benchmark tasks are cached on restart. The human data are not downloaded on every push.
