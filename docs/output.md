# Outputs

- `input_validation/`: normalized samplesheet and checksum manifest.
- `preprocessing/sample_objects/`: SeSAMe per-sample RDS and metrics; `preprocessing/cohort/`: compressed beta, M-value, detection matrices and analysis RDS.
- `qc/`: sample metrics, proposed/final exclusions, plots, and HTML QC report.
- `filtering/`: retained matrices, filtering waterfall, and every removed probe with reason.
- `differential_methylation/dmp/`: model definition plus complete/significant contrast tables.
- `report/`: cohort HTML and MultiQC custom-content table.
- `pipeline_info/`: Nextflow report, trace, timeline, DAG, and runtime provenance.
