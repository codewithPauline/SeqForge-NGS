process HAPLOTYPECALLER {
    tag { id }
    label 'calling'
    publishDir "${params.outdir}/variants/raw", mode: 'copy', pattern: '*.vcf.gz*'
    input:
    tuple val(id), path(bam), path(bai)
    path reference
    output:
    tuple val(id), path("${id}.raw.vcf.gz"), path("${id}.raw.vcf.gz.tbi"), emit: vcf
    path "${id}.gatk.versions.txt", emit: versions
    script:
    def heap = Math.max(1, task.memory.toGiga() - 1)
    """
    gatk --java-options '-Xmx${heap}g' HaplotypeCaller \
        -R '${reference}/reference.fa' -I '${bam}' \
        -L '${reference}/targets.bed' -O '${id}.raw.vcf.gz' \
        --sample-ploidy 2 --native-pair-hmm-threads ${task.cpus}
    gatk --version > '${id}.gatk.versions.txt' 2>&1
    """
}

process FILTER_VARIANTS {
    tag { id }
    publishDir "${params.outdir}/variants/filtered", mode: 'copy', pattern: '*.vcf.gz*'
    publishDir "${params.outdir}/qc/variants", mode: 'copy', pattern: '*.bcftools.stats'
    input:
    tuple val(id), path(vcf), path(tbi)
    path reference
    val minQual
    val minDepth
    output:
    tuple val(id), path("${id}.filtered.vcf.gz"), path("${id}.filtered.vcf.gz.tbi"), emit: filtered
    tuple val(id), path("${id}.pass.vcf.gz"), path("${id}.pass.vcf.gz.tbi"), emit: pass
    tuple val(id), path("${id}.bcftools.stats"), emit: stats
    path "${id}.bcftools.versions.txt", emit: versions
    script:
    """
    filter_variants.sh '${reference}/reference.fa' '${vcf}' '${id}' '${minQual}' '${minDepth}'
    bcftools --version > '${id}.bcftools.versions.txt'
    """
}
