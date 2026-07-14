# SeSAMe Nextflow Pipeline

This pipeline processes Illumina EPICv2 methylation array data (.idat files) using the [SeSAMe](https://www.bioconductor.org/packages/release/bioc/html/sesame.html) R package, wrapped in a Nextflow pipeline with Docker for environment control.

## Prerequisites
- Nextflow
- Docker

## Setup
1. **Build the Docker Image:**
   ```bash
   docker build -t sesame-pipeline:latest .
   ```

2. **Prepare Data:**
   Place your paired IDAT files (ending in `_Grn.idat` and `_Red.idat`) into a directory, e.g., `data/`.

## Running the Pipeline
Run the Nextflow pipeline, specifying the input pattern:

```bash
nextflow run main.nf --input "data/*_{Grn,Red}.idat" --outdir results
```

The pipeline will match the file pairs and process them using SeSAMe's `openSesame` function. The resulting beta values will be saved as CSV files in the specified `results` directory.
