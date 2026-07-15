# Feature parity assessment

Comparison baseline: [`nf-core/methylarray` development branch at `eb5fb7db5015069030bc6c0cf0d2f0985b93cefb`](https://github.com/nf-core/methylarray/tree/eb5fb7db5015069030bc6c0cf0d2f0985b93cefb), reviewed 2026-07-15.

| Feature | Current openSesame implementation | nf-core/methylarray behaviour | Proposed openSesame implementation | EPICv2 compatibility evidence | Implementation status | Known limitations |
|---|---|---|---|---|---|---|
| Input validation | File glob only | Samplesheet and schema | Validated CSV, normalized copy, checksums, deprecated glob converter | Platform-independent | Implemented | Local filesystem validation occurs in the workflow container context. |
| Platform verification | None | Assumes minfi-supported arrays | Per-IDAT SeSAMe read with recorded requested platform; cohort mixed-platform guard | SeSAMe supports multiple Infinium generations; platform is not inferred from probe count | Implemented | Exact platform classification depends on the installed SeSAMe manifest data. |
| Preprocessing | `openSesame()` beta CSV | minfi normalization | Explicit SeSAMe `readIDATpair()` + `prepSesame()` sequence, beta/M/detection outputs | SeSAMe release documentation describes preprocessing, QC and multi-generation array support | Implemented | The selected preparation code must be validated by the container's SeSAMe version. |
| Cohort matrices | Per-sample CSV | Cohort methylation objects | Compressed TSV matrices and `SummarizedExperiment` RDS | Matrix construction is platform-agnostic after verified preprocessing | Implemented | Annotation is conservative unless a vetted EPICv2 annotation is supplied. |
| Sample QC | None | minfi/ChAMP QC, plus incorrectly documented FastQC | Detection/intensity/missingness metrics, PCA, correlation, clustering, flags | SeSAMe supplies QC functions; no FastQC is appropriate for IDATs | Implemented | Sex/fingerprint calls are reported only when a validated installed method is available. |
| Probe filtering | None | Cross-reactive, SNP, sex and confounder filtering | Detection/missingness/user-list filtering with auditable reasons | EPICv2-specific masks must be versioned and supplied/verified | Partially implemented | Cross-reactive/SNP lists are intentionally not guessed; no phenotype-associated probe removal. |
| Cell composition | None | ChAMP reference adjustment | Disabled by default; validated user reference interface | No bundled reference is claimed EPICv2-valid | Deferred | Requires a verified EPICv2-compatible reference and overlap criterion. |
| Batch correction | None | Optional adjustment | QC assessment; model covariates; corrected values kept separate | Method is platform-independent but depends on valid design | Partially implemented | ComBat/removeBatchEffect execution is deferred pending a dedicated validated module. |
| DMP | None | ChAMP/minfi comparison | limma on filtered M-values with arbitrary formula and contrasts | limma is platform-independent when coordinates and inputs are valid | Implemented | Group replication and full-rank design are enforced. |
| DMR | None | DMR identification | DMRcate EPICv2 execution interface, BED/results outputs | DMRcate release documentation includes an EPICv2 vignette and suggests `EPICv2manifest` and `IlluminaHumanMethylationEPICv2anno.20a1.hg38` | Partially implemented | Requires those exact annotation packages in the execution image; otherwise stops clearly. |
| Blocks | None | Optional block analysis | Explicitly disabled/deferred module | No validated EPICv2 block method verified for this repository | Deferred | Blocks are not relabelled DMRs. |
| Aggregate report | None | MultiQC | Dedicated HTML cohort report and MultiQC custom-content input | Platform-independent | Implemented | MultiQC rendering requires MultiQC in the selected profile. |
| FastQC | None | Development README lists FastQC | Not implemented | IDATs are binary array intensity files, not sequencing reads | Not applicable | Deliberately excluded. |

This is a parity roadmap, not a claim of complete nf-core parity. Package and manifest versions used for every run are written to `pipeline_info/`.
