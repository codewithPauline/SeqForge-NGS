# Validation: software behavior and a regional human benchmark

The first HG002 calling/filtering benchmark is complete. See the [measured results](../benchmarks/evidence/hg002-grch38-chr20/README.md) and [reproduction instructions](benchmark.md). The protocol below records the broader acceptance criteria. Full FASTQ-to-VCF human benchmarking remains future work.

## What the current fixture establishes

The seeded generator creates a 12 kb artificial reference, two 3,000-pair FASTQ datasets, and a 10 kb target interval. SYNTH_A has one heterozygous SNV; SYNTH_B has a different homozygous alternate SNV. All sequence is synthetic, with fixed high base quality and no modeled sequencing errors. Gzip timestamps are fixed and generated files have a checksum manifest.

The integration checker verifies exact expected PASS variants and genotypes, sample names, BAM read groups and coordinate sorting, indexed regional access, mapped-read fraction, report counts, and report files. A second execution checks that every task is cached with `-resume`. The unit suite covers malformed reads, mate mismatches, contig errors, BED coordinates, and missing metrics.

These checks establish particular software behaviors on a simple input. They do not estimate performance on human reads, difficult regions, sequencing errors, or indels. Do not present the synthetic fixture as GIAB validation or turn its two successful calls into biological precision/recall claims.

## HG002 benchmark acceptance criteria

HG002 analyses require matching GIAB truth and confident regions. The completed regional benchmark starts from public alignments; future raw-read benchmarks will additionally evaluate upstream stages. Record:

| Resource | Required provenance |
|---|---|
| Sequencing reads | Provider, sample identity, library/platform, URLs/accessions, checksums, extraction or subsampling commands |
| Alignment reference | Exact assembly bundle, primary/alternate/decoy content, source URL, checksum, contig names |
| Truth VCF | HG002 identity, GIAB release, reference assembly, source URL, checksum |
| Confident-region BED | Same GIAB release and assembly as the truth VCF, URL and checksum |
| Evaluation region | Explicit autosomal coordinates; intersection with confident regions; any padding or excluded regions |
| Benchmark engine | Pinned hap.py or RTG vcfeval version, environment, exact invocation, and genotype/representation comparison settings |

Align reads to the full selected reference. Restrict evaluation to a documented autosomal interval and the confident regions. If extracting region-associated reads from an existing alignment, retain mate pairs and document that selection: it is a region-specific benchmark, not representative whole-genome validation. Boundary padding and the final scoring interval must be recorded separately.

Use a haplotype-aware benchmark tool to compare variants and genotypes. Simple coordinate overlap or `bcftools isec` alone is not sufficient for equivalent indel representations. Report SNP and INDEL results separately and retain the benchmark engine's original TP/FP/FN definitions, precision, recall, and F1 summaries.

Results require source manifests, commands, benchmark outputs, and scope limitations committed together. The current case study supplies those artifacts and uses preset thresholds. Future threshold selection and final evaluation should use different regions or samples to avoid reporting tuned performance as independent validation.

Starting resource index: [Genome in a Bottle data indexes](https://github.com/genome-in-a-bottle/giab_data_indexes). Selected URLs, object generations, sizes, and MD5 checksums are recorded in the [HG002 manifest](../benchmarks/hg002-grch38-chr20.json).
