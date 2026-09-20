# Changelog

## Unreleased

- Add a checksum-verified HG002/GRCh38 regional benchmark using GIAB v4.2.1 and RTG vcfeval 3.12.1, reusing the existing caller/filter modules.
- Publish measured SNP/INDEL metrics, source provenance, original ROC outputs, and false-positive/false-negative records for the fixed 1 Mb interval.
- Add benchmark verification, genotype/indel/region regression tests, and full benchmark integration/resume checks in CI.

- Add modular paired-end FASTQ validation, FastQC/fastp QC, BWA-MEM2 alignment, samtools duplicate marking, GATK HaplotypeCaller, and baseline bcftools filtering.
- Add MultiQC and SeqForge HTML/JSON reports, tool versions, reference provenance, and Nextflow execution reports.
- Add deterministic two-sample synthetic fixture, input-validation tests, expected-call integration checks, and resume verification.
- Add Docker build, exact Linux dependency lock, execution profiles, citation metadata, and methods/usage/benchmark documentation.

No tagged software release is published yet. Human accuracy metrics apply only to the documented regional calling/filtering benchmark.
