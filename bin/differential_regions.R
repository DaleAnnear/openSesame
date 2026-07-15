#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get <- function(key, default = NULL) {
  i <- match(key, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", key)
  args[i + 1]
}
required <- c("DMRcate", "DMRcatedata", "EPICv2manifest", "IlluminaHumanMethylationEPICv2anno.20a1.hg38", "limma")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("EPICv2 DMRcate requires missing package(s): ", paste(missing, collapse = ", "), ". Rebuild with the supplied Dockerfile.")

analysis <- readRDS(get("--analysis"))
assay <- function(x, name) if (inherits(x, "SummarizedExperiment")) SummarizedExperiment::assay(x, name) else x[[name]]
mvals <- assay(analysis, "mvalue")
if (!nrow(mvals)) stop("No probes remain after filtering; DMR analysis cannot proceed.")
metadata <- read.csv(get("--samplesheet"), check.names = FALSE, stringsAsFactors = FALSE)
metadata <- metadata[match(colnames(mvals), metadata$sample_id), , drop = FALSE]
if (anyNA(metadata$sample_id)) stop("Sample metadata does not match the filtered matrix.")

design <- model.matrix(as.formula(get("--design")), metadata)
if (qr(design)$rank < ncol(design)) stop("Design matrix is rank-deficient; DMRcate was not run.")
contrasts <- strsplit(get("--contrast"), ",", fixed = TRUE)[[1]]
contrasts <- trimws(contrasts[nzchar(trimws(contrasts))])
if (!length(contrasts)) stop("--contrast is required for DMR analysis.")
contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = design)

min_observed <- max(2L, ceiling(ncol(mvals) / 2))
mvals <- mvals[rowSums(is.finite(mvals)) >= min_observed, , drop = FALSE]
if (!nrow(mvals)) stop("No probes have methylation values in at least ", min_observed, " samples after filtering.")

outdir <- get("--outdir")
dir.create(file.path(outdir, "complete"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "significant"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
bandwidth <- as.numeric(get("--bandwidth"))
scaling <- as.numeric(get("--scaling"))
min_probes <- as.integer(get("--min-probes"))
fdr <- as.numeric(get("--fdr"))
min_delta <- as.numeric(get("--min-abs-delta-beta", "0"))
epic_filter <- get("--epicv2-filter", "mean")
if (!epic_filter %in% c("mean", "sensitivity", "precision", "random")) stop("--epicv2-filter must be mean, sensitivity, precision, or random.")
epic_remap <- tolower(get("--epicv2-remap", "true")) == "true"
build <- get("--genome-build", "hg38")
if (build != "hg38") stop("DMRcate EPICv2 support requires hg38 annotation; received ", build)

params <- data.frame(
  method = "DMRcate", arraytype = "EPICv2", genome_build = build,
  bandwidth = bandwidth, scaling = scaling, min_probes = min_probes,
  fdr = fdr, min_abs_delta_beta = min_delta, epicv2_filter = epic_filter,
  epicv2_remap = epic_remap, stringsAsFactors = FALSE
)
write.table(params, file.path(outdir, "dmr_parameters.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

summaries <- list()
for (contrast_name in colnames(contrast_matrix)) {
  safe <- gsub("[^A-Za-z0-9._-]", "_", contrast_name)
  cm <- contrast_matrix[, contrast_name, drop = FALSE]
  annotation <- DMRcate::cpg.annotate(
    datatype = "array", object = mvals, what = "M", arraytype = "EPICv2",
    epicv2Filter = epic_filter, epicv2Remap = epic_remap,
    analysis.type = "differential", design = design, contrasts = TRUE,
    cont.matrix = cm, coef = contrast_name, fdr = fdr
  )
  dmr_fit <- DMRcate::dmrcate(annotation, lambda = bandwidth, C = scaling)
  ranges <- DMRcate::extractRanges(dmr_fit, genome = build)
  complete <- as.data.frame(ranges)
  if (!nrow(complete)) {
    complete <- data.frame(seqnames = character(), start = integer(), end = integer(), no.cpgs = integer(), min_smoothed_fdr = numeric(), meandiff = numeric(), maxdiff = numeric())
  }
  complete$contrast <- contrast_name
  complete$representative_delta_beta <- if ("meandiff" %in% names(complete)) complete$meandiff else NA_real_
  complete$direction <- ifelse(complete$representative_delta_beta > 0, "hyper", ifelse(complete$representative_delta_beta < 0, "hypo", "neutral"))
  significant <- complete[complete$no.cpgs >= min_probes & complete$min_smoothed_fdr <= fdr & abs(complete$representative_delta_beta) >= min_delta, , drop = FALSE]
  write.table(complete, file.path(outdir, "complete", paste0(safe, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(significant, file.path(outdir, "significant", paste0(safe, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  bed <- data.frame(chrom = complete$seqnames, chromStart = pmax(0L, complete$start - 1L), chromEnd = complete$end, name = paste0("DMR_", seq_len(nrow(complete))), score = -log10(pmax(complete$min_smoothed_fdr, .Machine$double.xmin)), strand = ".")
  write.table(bed, file.path(outdir, "complete", paste0(safe, ".bed")), sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(list(annotation = annotation, dmrcate = dmr_fit, ranges = ranges), file.path(outdir, "complete", paste0(safe, ".rds")))
  if (nrow(significant)) {
    png(file.path(outdir, "plots", paste0(safe, "_top_dmr.png")), width = 1600, height = 1000, res = 160)
    groups <- if ("group" %in% names(metadata)) as.factor(metadata$group) else factor(rep("samples", nrow(metadata)))
    cols <- grDevices::hcl.colors(nlevels(groups), "Dark 3")[as.integer(groups)]
    tryCatch(DMRcate::DMR.plot(ranges = ranges, dmr = which(ranges$min_smoothed_fdr == min(ranges$min_smoothed_fdr))[1], CpGs = annotation, what = "Beta", arraytype = "EPICv2", phen.col = cols, genome = build), error = function(e) { plot.new(); text(.5, .5, paste("Top-DMR plot unavailable:", e$message)) })
    dev.off()
  }
  summaries[[contrast_name]] <- data.frame(contrast = contrast_name, tested_regions = nrow(complete), significant_regions = nrow(significant), stringsAsFactors = FALSE)
}
write.table(do.call(rbind, summaries), file.path(outdir, "dmr_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
