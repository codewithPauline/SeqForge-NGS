#!/usr/bin/env bash
# Shared by the Nextflow module and the filter boundary integration check.
set -euo pipefail
if [[ $# -ne 5 ]]; then
    echo 'Usage: filter_variants.sh reference.fa input.vcf sample min_qual min_depth' >&2
    exit 2
fi
reference=$1
input_vcf=$2
sample=$3
min_qual=$4
min_depth=$5
bcftools norm -f "$reference" -m -any "$input_vcf" -Ob -o normalized.bcf
bcftools filter --soft-filter LowQualDepth \
    --exclude "QUAL=\".\" || QUAL<${min_qual} || FMT/DP=\".\" || FMT/DP<${min_depth}" \
    normalized.bcf -Oz -o "${sample}.filtered.vcf.gz"
bcftools index --tbi "${sample}.filtered.vcf.gz"
bcftools view --apply-filters PASS "${sample}.filtered.vcf.gz" -Oz -o "${sample}.pass.vcf.gz"
bcftools index --tbi "${sample}.pass.vcf.gz"
bcftools stats "${sample}.pass.vcf.gz" > "${sample}.bcftools.stats"
