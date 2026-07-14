#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define default parameters
params.input = "data/*_{Grn,Red}.idat"
params.outdir = "results"

process processIdats {
    // Publish results to the specified output directory
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(idats)

    output:
    path "${sample_id}_betas.csv"

    script:
    // We expect the files to have the same prefix. The prefix is the part before _Grn.idat or _Red.idat.
    // It should be passed to the R script.
    // Because Nextflow stages files with their original names (if not renamed),
    // we can use sample_id as the prefix since it is the paired prefix.
    """
    process_idats.R --prefix ${sample_id} --output ${sample_id}_betas.csv
    """
}

workflow {
    // Input channel for IDAT pairs
    // The file pairing matches a common prefix and checks for both Grn and Red files
    Channel
        .fromFilePairs(params.input, size: 2)
        .set { idat_pairs_ch }

    processIdats(idat_pairs_ch)
}
