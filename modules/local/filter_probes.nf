process FILTER_PROBES {
    tag 'probe-filtering'
    label 'process_high'
    publishDir "${params.outdir}/filtering", mode: params.publish_mode
    input:
    path analysis
    path beta
    path mvalue
    path detection
    output:
    path 'filtering/analysis_filtered.rds', emit: analysis
    path 'filtering/beta_values.filtered.tsv.gz', emit: beta
    path 'filtering/m_values.filtered.tsv.gz', emit: mvalue
    path 'filtering/filtering_summary.tsv', emit: summary
    path 'filtering/probe_lists/removed_probes.tsv', emit: removed
    script:
    """filter_probes.R --analysis ${analysis} --outdir filtering --min-detection-rate ${params.min_detection_rate} --max-missing-rate ${params.max_missing_rate} --remove-sex ${params.remove_sex_chromosomes} --remove-noncpg ${params.remove_non_cpg} ${params.probe_exclusion_list ? "--exclusion-list ${params.probe_exclusion_list}" : ''}"""
}
