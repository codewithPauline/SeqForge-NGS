# Methods and design

This milestone is an engineering baseline for diploid, single-sample germline discovery. It does not implement the complete GATK Best Practices workflow.

1. Validate every samplesheet row before scheduling reference indexing or sample tasks. Stream both FASTQ mates to check record structure, pair identifiers, and counts.
2. Measure raw-read QC with FastQC; trim and filter paired reads with fastp using paired-end adapter detection. Preserve the fastp before/after metrics.
3. Validate FASTA/BED contigs and bounds. Hash the reference and effective BED, and create BWA-MEM2, FASTA, and sequence-dictionary indexes.
4. Align paired reads with BWA-MEM2 and explicit Illumina read groups. Pipe SAM output into name sorting, run `samtools fixmate -m`, coordinate sort, then `samtools markdup`. Duplicates are **marked and retained**, not removed. Index and check the BAM.
5. Capture `samtools stats`, `flagstat`, and duplicate-marking metrics. Call variants with GATK HaplotypeCaller in ordinary VCF mode, diploid, over the effective intervals.
6. Normalize against the same FASTA and split multiallelic records with bcftools. Mark low/missing QUAL or sample depth; preserve both labeled and PASS-only files. Compute statistics from PASS calls.
7. Aggregate QC with MultiQC and generate a separate HTML/JSON summary. Capture executable version output, reference and interval checksums, and Nextflow execution records.

The reference is shared through a Nextflow value channel, so every sample receives the same indexed reference. Per-sample joins use the validated sample ID. Input basename restrictions keep filenames safe in generated shell commands. The CI fixture uses two samples with distinct variant/genotype expectations to expose sample mixing or accidental reference-channel consumption.

## Deliberate limitations

- No BQSR or known-sites resource management yet.
- No cohort gVCF/joint genotyping, VQSR, GATK SNP/indel-specific hard filtering, or optimized exome filtering.
- No non-diploid or sex-aware calling, structural variant calling, CNV calling, or somatic analysis.
- No functional annotation yet. VEP requires a versioned cache and assembly-matched resources.
- No read-depth coverage gate, contamination estimate, sample fingerprinting, or clinical interpretation.
- No simulated INDEL test in this first fixture: the caller supports small indels, but this milestone's expected-call integration assertions exercise SNVs only.
- A separate HG002 regional benchmark now evaluates calling/filtering from public alignments. Full upstream human validation and measured cluster-scale performance remain pending.

## Primary documentation

- [Nextflow process and dataflow semantics](https://www.nextflow.io/docs/stable/process.html)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [fastp](https://github.com/OpenGene/fastp)
- [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)
- [samtools markdup prerequisites and behavior](https://www.htslib.org/doc/samtools-markdup.html)
- [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
- [bcftools normalization, expressions, and filtering](https://samtools.github.io/bcftools/bcftools.html)
- [MultiQC](https://docs.seqera.io/multiqc/)

Tool versions are specified in `environment.yml` and exact Linux builds in `environment-linux-64.lock`. Cite the corresponding tools as well as SeqForge-NGS in scientific work.
