process COHORT_REPORT {
    tag 'report'
    label 'process_low'
    publishDir "${params.outdir}/report", mode: params.publish_mode
    input:
    path manifest
    path qc
    path exclusions
    path filtering_summary
    output:
    path 'report/cohort_report.html', emit: report
    path 'report/multiqc_opensesame_general_stats.tsv', emit: multiqc
    script:
    """cohort_report.R --manifest ${manifest} --qc ${qc} --exclusions ${exclusions} --filtering ${filtering_summary} --outdir report"""
}
