process PREPROCESS_IDAT {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/preprocessing/sample_objects", mode: params.publish_mode
    input: tuple val(sample_id), path(red), path(green), val(clean_id)
    output:
    tuple val(sample_id), path("${clean_id}.rds"), emit: objects
    path "${clean_id}.metrics.tsv", emit: metrics
    path "${clean_id}.versions.tsv", emit: versions
    script:
    """preprocess_idat.R --sample-id '${sample_id}' --red ${red} --green ${green} --prep-code '${params.sesame_prep_code}' --output ${clean_id}.rds --metrics ${clean_id}.metrics.tsv --versions ${clean_id}.versions.tsv"""
}
