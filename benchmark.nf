nextflow.enable.dsl = 2
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import java.security.MessageDigest

include { HAPLOTYPECALLER; FILTER_VARIANTS } from './modules/variants'
include { VERIFY_BENCHMARK; RTG_FORMAT; RTG_COMPARE; BENCHMARK_REPORT } from './modules/benchmark'

workflow {
    if (!params.benchmark_data) error 'Required: --benchmark_data path/to/prepared/bundle. See docs/benchmark.md.'
    def bundle = file(params.benchmark_data, checkIfExists: true)
    def prepared = bundle.resolve('prepared.json')
    def metadata = new JsonSlurper().parse(prepared.toFile())
    def manifest = metadata.manifest
    if (!(manifest.sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]{0,79}/)) error 'Invalid benchmark sample ID'
    if (!(manifest.min_depth.toString() ==~ /[1-9][0-9]*/)) error 'Invalid benchmark minimum depth'
    if (!(manifest.min_qual.toString() ==~ /[0-9]+(\.[0-9]+)?/)) error 'Invalid benchmark minimum QUAL'
    if (!(params.max_cpus.toString() ==~ /[1-9][0-9]*/)) error 'Invalid --max_cpus'
    def sourcePaths = ['benchmark.nf', 'nextflow.config', 'modules/variants.nf', 'modules/benchmark.nf',
        'bin/filter_variants.sh', 'bin/verify_benchmark.py', 'bin/benchmark_io.py', 'bin/benchmark_report.py',
        'scripts/prepare_benchmark.py', 'environment-linux-64.lock']
    def source = sourcePaths.collectEntries { relative ->
        [(relative): MessageDigest.getInstance('SHA-256').digest(projectDir.resolve(relative).bytes).encodeHex().toString()]
    }
    VERIFY_BENCHMARK(Channel.value(bundle))
    reference = VERIFY_BENCHMARK.out.verified.map { bundle.resolve('reference') }
    reads = VERIFY_BENCHMARK.out.verified.map { tuple(manifest.sample, bundle.resolve('reads.bam'), bundle.resolve('reads.bam.bai')) }
    HAPLOTYPECALLER(reads, reference)
    FILTER_VARIANTS(HAPLOTYPECALLER.out.vcf, reference, manifest.min_qual, manifest.min_depth)
    RTG_FORMAT(reference)
    RTG_COMPARE(FILTER_VARIANTS.out.pass, RTG_FORMAT.out.sdf,
        Channel.value(bundle.resolve('truth.vcf.gz')), Channel.value(bundle.resolve('truth.vcf.gz.tbi')),
        Channel.value(bundle.resolve('evaluation.bed')))
    versions = HAPLOTYPECALLER.out.versions.mix(FILTER_VARIANTS.out.versions)
        .mix(RTG_FORMAT.out.versions).mix(RTG_COMPARE.out.versions).collect(sort: true)
    BENCHMARK_REPORT(RTG_COMPARE.out.results, Channel.value(prepared),
        FILTER_VARIANTS.out.pass.map { id, vcf, tbi -> vcf }, versions, JsonOutput.toJson(source))
}
