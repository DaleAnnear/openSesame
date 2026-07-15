process DIFFERENTIAL_REGIONS {
    tag 'dmr'
    label 'process_high'
    publishDir "${params.outdir}/differential_methylation/dmr", mode: params.publish_mode
    input:
    path complete
    path analysis
    path samplesheet
    output: path 'dmr', emit: results
    script:
    """differential_regions.R --dmp-dir ${complete} --analysis ${analysis} --outdir dmr --bandwidth ${params.dmr_bandwidth} --scaling ${params.dmr_scaling} --min-probes ${params.dmr_min_probes}"""
}
