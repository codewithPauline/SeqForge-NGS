process ALIGN_MARKDUP {
    tag { id }
    label 'align'
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: '*.{bam,bai}'
    publishDir "${params.outdir}/qc/alignment", mode: 'copy', pattern: '*.markdup.txt'
    input:
    tuple val(id), path(r1), path(r2)
    path reference
    output:
    tuple val(id), path("${id}.bam"), path("${id}.bam.bai"), emit: bam
    tuple val(id), path("${id}.markdup.txt"), emit: metrics
    path "${id}.alignment.versions.txt", emit: versions
    script:
    // Reserve one core for the simultaneous sorter. -@ counts EXTRA samtools threads.
    def alignThreads = Math.max(1, task.cpus - 1)
    """
    bwa-mem2 mem -t ${alignThreads} \
        -R '@RG\\tID:${id}\\tSM:${id}\\tPL:ILLUMINA\\tLB:${id}\\tPU:${id}' \
        '${reference}/reference.fa' '${r1}' '${r2}' \
        | samtools sort -n -@ 0 -m 256M -o namesort.bam -
    samtools fixmate -m -@ 0 namesort.bam fixmate.bam
    samtools sort -@ 0 -m 512M -o coordsort.bam fixmate.bam
    samtools markdup -@ 0 -s -f '${id}.markdup.txt' coordsort.bam '${id}.bam'
    samtools index '${id}.bam'
    samtools quickcheck -v '${id}.bam'
    rm namesort.bam fixmate.bam coordsort.bam
    bwa-mem2 version > '${id}.alignment.versions.txt' 2>&1
    samtools --version >> '${id}.alignment.versions.txt'
    """
}

process BAM_QC {
    tag { id }
    publishDir "${params.outdir}/qc/alignment", mode: 'copy', pattern: '*.{stats,flagstat.txt}'
    input:
    tuple val(id), path(bam), path(bai)
    output:
    tuple val(id), path("${id}.samtools.stats"), path("${id}.flagstat.txt"), emit: qc
    path "${id}.bamqc.versions.txt", emit: versions
    script:
    """
    samtools stats '${bam}' > '${id}.samtools.stats'
    samtools flagstat '${bam}' > '${id}.flagstat.txt'
    samtools --version > '${id}.bamqc.versions.txt'
    """
}
