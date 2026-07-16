#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VALIDATE_SAMPLESHEET } from './modules/local/validate_samplesheet'
include { CREATE_LEGACY_SAMPLESHEET } from './modules/local/create_legacy_samplesheet'
include { PREPROCESS_IDAT } from './modules/local/preprocess_idat'
include { BUILD_COHORT } from './modules/local/build_cohort'
include { SAMPLE_QC } from './modules/local/sample_qc'
include { FILTER_PROBES } from './modules/local/filter_probes'
include { DIFFERENTIAL_METHYLATION } from './modules/local/differential_methylation'
include { DIFFERENTIAL_REGIONS } from './modules/local/differential_regions'
include { COHORT_REPORT } from './modules/local/cohort_report'

workflow {
    if( !params.input && !params.input_glob ) error "Provide --input <samplesheet.csv> (preferred) or deprecated --input_glob."
    if( params.input && params.input_glob ) error "Use only one of --input and --input_glob."
    samplesheet_ch = params.input ? Channel.fromPath(params.input, checkIfExists: true) : null
    if( params.input_glob ) {
        log.warn "--input_glob is deprecated; use a samplesheet with sample_id,idat_red,idat_green."
        CREATE_LEGACY_SAMPLESHEET(Channel.fromPath(params.input_glob, checkIfExists: true).collect())
        samplesheet_ch = CREATE_LEGACY_SAMPLESHEET.out.samplesheet
    }
    VALIDATE_SAMPLESHEET(samplesheet_ch)
    validated_ch = VALIDATE_SAMPLESHEET.out.samplesheet
    idat_ch = validated_ch.splitCsv(header: true).map { row ->
        def fields = row.collectEntries { key, entry -> [(key.toString().trim().replaceAll(/^"|"$/, '')): entry?.toString()?.trim()?.replaceAll(/^"|"$/, '')] }
        def sampleId = fields['sample_id']; def redIdat = fields['idat_red']; def greenIdat = fields['idat_green']
        if( !sampleId || !redIdat || !greenIdat ) error "Validated samplesheet row is missing sample_id, idat_red, or idat_green: ${fields}"
        def clean = sampleId.replaceAll(/[^A-Za-z0-9._-]/, '_')
        tuple(sampleId, file(redIdat), file(greenIdat), clean)
    }
    PREPROCESS_IDAT(idat_ch)
    sample_objects_ch = PREPROCESS_IDAT.out.objects.map { sample_id, object_file -> object_file }.collect()
    BUILD_COHORT(sample_objects_ch, validated_ch)
    SAMPLE_QC(BUILD_COHORT.out.analysis, BUILD_COHORT.out.beta, BUILD_COHORT.out.mvalue, BUILD_COHORT.out.detection, validated_ch)
    FILTER_PROBES(BUILD_COHORT.out.analysis, BUILD_COHORT.out.beta, BUILD_COHORT.out.mvalue, BUILD_COHORT.out.detection)
    run_dmps = params.find_dmps.toString().toBoolean()
    run_dmrs = params.find_dmrs.toString().toBoolean()
    if( run_dmps ) {
        DIFFERENTIAL_METHYLATION(FILTER_PROBES.out.analysis, FILTER_PROBES.out.beta, FILTER_PROBES.out.mvalue, validated_ch)
        if( run_dmrs ) DIFFERENTIAL_REGIONS(DIFFERENTIAL_METHYLATION.out.complete, FILTER_PROBES.out.analysis, validated_ch)
    } else if( run_dmrs ) error "--find_dmrs requires --find_dmps true because regions reuse the validated DMP model."
    COHORT_REPORT(VALIDATE_SAMPLESHEET.out.manifest, SAMPLE_QC.out.qc, SAMPLE_QC.out.exclusions, FILTER_PROBES.out.summary)
}
