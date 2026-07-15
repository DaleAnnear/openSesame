process CREATE_LEGACY_SAMPLESHEET {
    tag 'legacy-idat-glob'
    label 'process_low'
    input: path idats
    output: path 'legacy.samplesheet.csv', emit: samplesheet
    script:
    """create_legacy_samplesheet.R --output legacy.samplesheet.csv ${idats}"""
}
