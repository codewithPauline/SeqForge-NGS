process VERIFY_BENCHMARK {
    input:
    path bundle, stageAs: 'bundle'
    output:
    path 'verification.json', emit: verified
    script:
    """
    verify_benchmark.py bundle
    """
}

process RTG_FORMAT {
    memory '4 GB'
    input:
    path reference
    output:
    path 'reference.sdf', emit: sdf
    path 'rtg_format.versions.txt', emit: versions
    script:
    """
    RTG_MEM=3g rtg format -o reference.sdf '${reference}/reference.fa'
    rtg version > rtg_format.versions.txt
    """
}

process RTG_COMPARE {
    label 'calling'
    publishDir "${params.outdir}/benchmark", mode: 'copy', pattern: 'comparison'
    input:
    tuple val(id), path(vcf), path(tbi)
    path sdf
    path truth
    path truth_index
    path evaluation
    output:
    path 'comparison', emit: results
    path 'rtg_compare.versions.txt', emit: versions
    script:
    def heap = Math.max(1, task.memory.toGiga() - 1)
    """
    RTG_MEM=${heap}g rtg vcfeval --baseline '${truth}' --calls '${vcf}' \
        --template '${sdf}' --evaluation-regions '${evaluation}' --output comparison \
        --sample '${id}' --sample-ploidy 2 --vcf-score-field QUAL \
        --roc-subset snp,indel --threads ${task.cpus}
    rtg version > rtg_compare.versions.txt
    """
}

process BENCHMARK_REPORT {
    publishDir "${params.outdir}/benchmark", mode: 'copy', pattern: 'benchmark_summary.*'
    publishDir "${params.outdir}/benchmark", mode: 'copy', pattern: 'prepared.json'
    publishDir "${params.outdir}/benchmark", mode: 'copy', pattern: 'source_snapshot.json'
    input:
    path comparison
    path prepared, stageAs: 'prepared.json'
    path query
    path versions, stageAs: 'versions/*'
    val source_json
    output:
    path 'benchmark_summary.json', emit: json
    path 'benchmark_summary.html', emit: html
    path 'prepared.json', emit: provenance
    path 'source_snapshot.json', emit: source
    script:
    """
    echo '${source_json}' > source_snapshot.json
    benchmark_report.py --comparison '${comparison}' --prepared prepared.json \
        --query '${query}' --versions versions --source source_snapshot.json
    """
}
