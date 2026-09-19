process MULTIQC {
    publishDir "${params.outdir}/reports", mode: 'copy', pattern: 'multiqc*'
    input:
    path qc_files, stageAs: 'qc/*'
    output:
    path 'multiqc_report.html', emit: html
    path 'multiqc_report_data', emit: data
    path 'multiqc.versions.txt', emit: versions
    script:
    """
    multiqc qc --force --filename multiqc_report.html --title 'SeqForge-NGS QC'
    multiqc --version > multiqc.versions.txt
    """
}

process RUN_SUMMARY {
    publishDir "${params.outdir}/reports", mode: 'copy', pattern: 'seqforge*'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', pattern: 'software_versions.txt'
    input:
    path metrics, stageAs: 'metrics/*'
    path versions, stageAs: 'versions/*'
    path metadata
    val pipeline_version
    val nextflow_version
    val minQual
    val minDepth
    output:
    path 'seqforge_summary.json', emit: json
    path 'seqforge_summary.html', emit: html
    path 'software_versions.txt', emit: software
    script:
    """
    make_report.py --metrics metrics --versions versions --reference '${metadata}' \
        --pipeline-version '${pipeline_version}' --nextflow-version '${nextflow_version}' \
        --min-qual '${minQual}' --min-depth '${minDepth}'
    """
}
