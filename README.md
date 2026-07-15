# openSesame

`openSesame` is a modular Nextflow DSL2 workflow for Illumina methylation IDATs, using SeSAMe as its primary preprocessing engine. Its supported target is MethylationEPIC v2.0 (EPICv2); support is deliberately conservative and is recorded in [the feature-parity matrix](docs/feature_parity.md).

## Quick start

Build the pinned image, then run with a validated samplesheet:

```bash
docker build -t ghcr.io/daleannear/opensesame:0.2.0 .
nextflow run . --input samplesheet.csv --outdir results -profile docker
```

For an association analysis:

```bash
nextflow run . --input samplesheet.csv --outdir results \
  --design '~ 0 + group + age + sex + batch' \
  --contrast 'groupDisease-groupControl' --find_dmps true -profile docker
```

Apptainer and generic SLURM profiles are also supplied: `-profile apptainer` and `-profile slurm,apptainer`. Site resource overrides belong in a separate `-c` config. JSON/YAML parameter files are supported through Nextflow `-params-file`.

The input CSV must contain `sample_id,idat_red,idat_green`; arbitrary additional metadata are retained for modelling. Both `.idat` and `.idat.gz` are accepted. `--input_glob` remains temporarily available but is deprecated.

## Outputs and safeguards

The pipeline writes validated inputs and provenance under `pipeline_info/` and `input_validation/`, unfiltered cohort matrices under `preprocessing/`, QC under `qc/`, auditable filtering under `filtering/`, and DMP/DMR results under `differential_methylation/`. Beta values are for interpretation; M-values are used by limma. Flagged samples are not removed unless `--auto_exclude_samples true` is set. Unfiltered and uncorrected matrices are retained.

See [usage](docs/usage.md), [outputs](docs/output.md), [methods](docs/methods.md), and [troubleshooting](docs/troubleshooting.md). Run `nextflow config .`, `nf-test test`, and the test profile in an environment with Nextflow, nf-test, Docker, and the documented public EPICv2 fixture installed.
