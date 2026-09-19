nextflow.enable.dsl = 2

include { PREPARE_REFERENCE } from './modules/reference'
include { CHECK_FASTQ; FASTQC; FASTP } from './modules/reads'
include { ALIGN_MARKDUP; BAM_QC } from './modules/alignment'
include { HAPLOTYPECALLER; FILTER_VARIANTS } from './modules/variants'
include { MULTIQC; RUN_SUMMARY } from './modules/reporting'

workflow {
    if (params.help) {
        log.info '''
SeqForge-NGS — reproducible paired-end germline variant discovery

  nextflow run main.nf -profile docker \
    --input samplesheet.csv --fasta reference.fa --outdir results

Required: --input CSV (sample,fastq_1,fastq_2), --fasta uncompressed FASTA
Optional: --intervals BED, --max_cpus 8, --min_qual 30, --min_depth 10
          --container seqforge-ngs:0.1.0-dev
Profiles: docker, singularity, slurm, test (combine: -profile test,docker)
Test data: python scripts/make_test_data.py
Details: README.md and docs/usage.md
'''
        return
    }
    if (!params.input || !params.fasta)
        error 'Required: --input samplesheet.csv --fasta reference.fa. See --help.'
    if (!(params.max_cpus.toString() ==~ /[1-9][0-9]*/))
        error '--max_cpus must be a positive integer'
    if (!(params.min_depth.toString() ==~ /[1-9][0-9]*/))
        error '--min_depth must be a positive integer'
    if (!(params.min_qual.toString() ==~ /[0-9]+(\.[0-9]+)?/))
        error '--min_qual must be a nonnegative number'

    def sheet = Inputs.localFile(params.input, launchDir, 'Samplesheet')
    def fasta = Inputs.localFile(params.fasta, launchDir, 'Reference')
    if (fasta.toString().endsWith('.gz')) error 'Reference must be uncompressed FASTA'
    def bed = params.intervals ? Inputs.localFile(params.intervals, launchDir, 'Intervals') : []

    // Validate the entire sheet before emitting any sample to expensive processes.
    validated_samples = Channel.fromPath(sheet).splitCsv(header: true, strip: true)
        .toList().map { rows -> Inputs.samples(rows, sheet.parent) }
    samples = validated_samples.flatMap { it }
    PREPARE_REFERENCE(validated_samples.map { fasta }, Channel.value(bed))
    CHECK_FASTQ(samples)
    FASTQC(CHECK_FASTQ.out.reads)
    FASTP(CHECK_FASTQ.out.reads)
    ALIGN_MARKDUP(FASTP.out.reads, PREPARE_REFERENCE.out.bundle)
    BAM_QC(ALIGN_MARKDUP.out.bam)
    HAPLOTYPECALLER(ALIGN_MARKDUP.out.bam, PREPARE_REFERENCE.out.bundle)
    FILTER_VARIANTS(HAPLOTYPECALLER.out.vcf, PREPARE_REFERENCE.out.bundle,
        params.min_qual, params.min_depth)

    qc_files = FASTQC.out.zip.map { id, files -> files }
        .mix(FASTP.out.json.map { id, f -> f })
        .mix(BAM_QC.out.qc.map { id, stats, flagstat -> [stats, flagstat] })
        .mix(ALIGN_MARKDUP.out.metrics.map { id, f -> f })
        .mix(FILTER_VARIANTS.out.stats.map { id, f -> f })
        .collect()
    MULTIQC(qc_files)
    summaries = FASTP.out.json.join(BAM_QC.out.qc).join(FILTER_VARIANTS.out.stats)
        .map { id, fastp, stats, flagstat, variants -> [fastp, stats, flagstat, variants] }
        .collect()
    versions = PREPARE_REFERENCE.out.versions.mix(CHECK_FASTQ.out.versions)
        .mix(FASTQC.out.versions).mix(FASTP.out.versions)
        .mix(ALIGN_MARKDUP.out.versions).mix(BAM_QC.out.versions)
        .mix(HAPLOTYPECALLER.out.versions).mix(FILTER_VARIANTS.out.versions)
        .mix(MULTIQC.out.versions).collect()
    RUN_SUMMARY(summaries, versions, PREPARE_REFERENCE.out.metadata,
        workflow.manifest.version, workflow.nextflow.version.toString(),
        params.min_qual, params.min_depth)
}
