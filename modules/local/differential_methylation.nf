process DIFFERENTIAL_METHYLATION {
    tag 'dmp'
    label 'process_high'
    publishDir "${params.outdir}/differential_methylation/dmp", mode: params.publish_mode
    input:
    path analysis
    path beta
    path mvalue
    path samplesheet
    output:
    path 'dmp/complete', emit: complete
    path 'dmp/significant', emit: significant
    path 'dmp/model_definitions.tsv', emit: model
    path 'dmp/dmp_testability_summary.tsv', emit: testability
    script:
    """differential_methylation.R --implementation-version 4 --analysis ${analysis} --samplesheet ${samplesheet} --outdir dmp --design '${params.design}' --contrast '${params.contrast}' --fdr ${params.fdr} --min-abs-delta-beta ${params.min_abs_delta_beta}"""
}
