# EPICv2 DMR analysis

When `--find_dmrs true`, openSesame runs DMRcate on filtered M-values using the same `--design` and comma-separated `--contrast` values as DMP analysis. This is an EPICv2-specific implementation: DMRcate uses the EPICv2 manifest to make the required one-to-one probe-to-CpG representation before smoothing. The default `--dmr_epicv2_filter mean` averages duplicate-target probes; `sensitivity` and `precision` are alternative documented selection strategies, while `random` is not reproducible and should generally be avoided.

The process uses hg38 only, requires `DMRcate`, `DMRcatedata`, `EPICv2manifest`, and `IlluminaHumanMethylationEPICv2anno.20a1.hg38`, and stops rather than using EPICv1 coordinates.

For each contrast, `differential_methylation/dmr/complete/` contains all regions as TSV, BED, and RDS; `significant/` contains regions meeting the smoothed FDR, minimum CpGs, and `--min_abs_delta_beta` criteria. `representative_delta_beta` is DMRcate's regional mean beta difference (`meandiff`); `maxdiff` is retained as the maximum difference. Top DMR plots, a summary, and all DMR parameters are also written.
