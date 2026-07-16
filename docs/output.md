# Outputs

- `input_validation/`: normalized samplesheet and checksum manifest.
- `preprocessing/sample_objects/`: SeSAMe per-sample RDS and metrics; `preprocessing/cohort/`: compressed beta, M-value, detection matrices and analysis RDS.
- `qc/`: sample metrics, proposed/final exclusions, plots, and HTML QC report.
- `filtering/`: retained matrices, filtering waterfall, and every removed probe with reason.
- `differential_methylation/dmp/`: model definition, testability summary, and complete/significant contrast tables. For a two-group contrast such as `groupCase-groupControl`, each complete DMP table includes `delta_beta_adjusted` (the covariate-adjusted, unconstrained beta-scale model coefficient) and `delta_beta_observed` (the observed mean beta in Case minus Control, bounded to -1 through 1). The significant table applies `--min_abs_delta_beta` to `delta_beta_observed`.
- `report/`: cohort HTML and MultiQC custom-content table.
- `pipeline_info/`: Nextflow report, trace, timeline, DAG, and runtime provenance.
