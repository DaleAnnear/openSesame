process SAMPLE_QC {
    tag 'cohort-qc'
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: params.publish_mode
    input:
    path analysis
    path beta
    path mvalue
    path detection
    path samplesheet
    output:
    path 'qc/sample_qc.tsv', emit: qc
    path 'qc/proposed_exclusions.tsv', emit: exclusions
    path 'qc/qc_report.html', emit: report
    script:
    """sample_qc.R --analysis ${analysis} --samplesheet ${samplesheet} --outdir qc --min-detection-rate ${params.min_detection_rate} --auto-exclude ${params.auto_exclude_samples}"""
}
