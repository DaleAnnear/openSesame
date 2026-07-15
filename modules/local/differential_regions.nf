process DIFFERENTIAL_REGIONS {
    tag 'dmrcate-epicv2'
    label 'process_high'
    publishDir "${params.outdir}/differential_methylation/dmr", mode: params.publish_mode
    input:
    path complete
    path analysis
    path samplesheet
    output:
    path 'dmr/complete', emit: complete
    path 'dmr/significant', emit: significant
    path 'dmr/dmr_summary.tsv', emit: summary
    path 'dmr/dmr_parameters.tsv', emit: parameters
    path 'dmr/plots', emit: plots
    script:
    """
    differential_regions.R --dmp-dir ${complete} --analysis ${analysis} --samplesheet ${samplesheet} --outdir dmr --design '${params.design}' --contrast '${params.contrast}' --bandwidth ${params.dmr_bandwidth} --scaling ${params.dmr_scaling} --min-probes ${params.dmr_min_probes} --fdr ${params.dmr_fdr} --min-abs-delta-beta ${params.min_abs_delta_beta} --epicv2-filter '${params.dmr_epicv2_filter}' --epicv2-remap ${params.dmr_epicv2_remap} --genome-build ${params.genome_build}
    """
}
