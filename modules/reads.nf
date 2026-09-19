process CHECK_FASTQ {
    tag { id }
    publishDir "${params.outdir}/qc/input", mode: 'copy', pattern: '*.input.json'
    input:
    tuple val(id), path(r1, stageAs: 'r1/*'), path(r2, stageAs: 'r2/*')
    output:
    tuple val(id), path(r1), path(r2), emit: reads
    path "${id}.input.json", emit: validation
    path "${id}.input.versions.txt", emit: versions
    script:
    """
    validate_fastq.py --read1 '${r1}' --read2 '${r2}' --sample '${id}' \
        --output '${id}.input.json'
    python --version > '${id}.input.versions.txt'
    """
}

process FASTQC {
    tag { id }
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy', pattern: '*_fastqc.*'
    input:
    tuple val(id), path(r1, stageAs: 'r1/*'), path(r2, stageAs: 'r2/*')
    output:
    tuple val(id), path('*_fastqc.zip'), emit: zip
    tuple val(id), path('*_fastqc.html'), emit: html
    path "${id}.fastqc.versions.txt", emit: versions
    script:
    def suffix1 = r1.name.endsWith('.gz') ? 'fastq.gz' : 'fastq'
    def suffix2 = r2.name.endsWith('.gz') ? 'fastq.gz' : 'fastq'
    """
    ln -s '${r1}' '${id}_R1.${suffix1}'
    ln -s '${r2}' '${id}_R2.${suffix2}'
    fastqc --threads 1 --outdir . '${id}_R1.${suffix1}' '${id}_R2.${suffix2}'
    fastqc --version > '${id}.fastqc.versions.txt'
    """
}

process FASTP {
    tag { id }
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', pattern: '*.{json,html}'
    input:
    tuple val(id), path(r1, stageAs: 'r1/*'), path(r2, stageAs: 'r2/*')
    output:
    tuple val(id), path("${id}_R1.trim.fastq.gz"), path("${id}_R2.trim.fastq.gz"), emit: reads
    tuple val(id), path("${id}.fastp.json"), emit: json
    tuple val(id), path("${id}.fastp.html"), emit: html
    path "${id}.fastp.versions.txt", emit: versions
    script:
    """
    fastp --in1 '${r1}' --in2 '${r2}' \
        --out1 '${id}_R1.trim.fastq.gz' --out2 '${id}_R2.trim.fastq.gz' \
        --detect_adapter_for_pe --thread ${task.cpus} \
        --json '${id}.fastp.json' --html '${id}.fastp.html'
    fastp --version > '${id}.fastp.versions.txt' 2>&1
    """
}
