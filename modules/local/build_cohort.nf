process BUILD_COHORT {
    tag 'cohort'
    label 'process_high'
    publishDir "${params.outdir}/preprocessing", mode: params.publish_mode
    input:
    path objects
    path samplesheet
    output:
    path 'cohort/analysis_unfiltered.rds', emit: analysis
    path 'cohort/beta_values.tsv.gz', emit: beta
    path 'cohort/m_values.tsv.gz', emit: mvalue
    path 'cohort/detection_metrics.tsv.gz', emit: detection
    path 'cohort/probe_annotation.tsv.gz', emit: annotation
    path 'cohort/sample_metadata.tsv', emit: metadata
    script:
    """build_cohort.R --samplesheet ${samplesheet} --outdir cohort ${objects}"""
}
