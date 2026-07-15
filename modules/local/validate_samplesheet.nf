process VALIDATE_SAMPLESHEET {
    tag 'samplesheet'
    label 'process_low'
    publishDir "${params.outdir}/input_validation", mode: params.publish_mode
    input: path samplesheet
    output:
    path 'samplesheet.valid.csv', emit: samplesheet
    path 'run_manifest.json', emit: manifest
    script:
    """validate_samplesheet.R --input ${samplesheet} --output samplesheet.valid.csv --manifest run_manifest.json"""
}
