process PREPARE_REFERENCE {
    label 'index'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'reference.json'
    input:
    path fasta, stageAs: 'input/reference.fa'
    path intervals, stageAs: 'input/targets.bed'
    output:
    path 'reference', emit: bundle
    path 'reference.json', emit: metadata
    path 'reference.versions.txt', emit: versions
    script:
    def bedArg = intervals ? '--intervals input/targets.bed' : ''
    """
    mkdir reference
    validate_reference.py --fasta input/reference.fa ${bedArg} \
        --bed-out reference/targets.bed --metadata reference.json
    cp input/reference.fa reference/reference.fa
    samtools faidx reference/reference.fa
    samtools dict -o reference/reference.dict reference/reference.fa
    bwa-mem2 index reference/reference.fa
    bwa-mem2 version > reference.versions.txt 2>&1
    samtools --version >> reference.versions.txt
    """
}
