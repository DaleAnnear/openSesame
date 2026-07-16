#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get <- function(k) { i <- match(k, args); if (is.na(i) || i == length(args)) stop("Missing ", k); args[i + 1] }
if (!requireNamespace("limma", quietly = TRUE)) stop("limma is required for DMP analysis")

a <- readRDS(get("--analysis"))
ss <- read.csv(get("--samplesheet"), check.names = FALSE, stringsAsFactors = FALSE)
out <- get("--outdir")
dir.create(file.path(out, "complete"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "significant"), recursive = TRUE, showWarnings = FALSE)
assay <- function(name) if (inherits(a, "SummarizedExperiment")) SummarizedExperiment::assay(a, name) else a[[name]]
m <- assay("mvalue")
b <- assay("beta")
if (!nrow(m)) stop("No probes remain after filtering. Inspect filtering/filtering_summary.tsv and filtering/probe_lists/removed_probes.tsv; for small cohorts, consider scientifically justified relaxation of --min_detection_rate and/or --max_missing_rate.")
if (ncol(m) < 2) stop("DMP analysis requires at least two samples after QC.")
if (!identical(rownames(m), rownames(b)) || !identical(colnames(m), colnames(b))) stop("M-value and beta matrices are not aligned.")

ss <- ss[match(colnames(m), ss$sample_id), , drop = FALSE]
if (anyNA(ss$sample_id)) stop("Samplesheet does not match the analysis matrix.")
design <- model.matrix(as.formula(get("--design")), ss)
if (qr(design)$rank < ncol(design)) stop("Design matrix is rank-deficient; revise design/covariates")
contrast <- get("--contrast")
if (contrast == "" || is.na(contrast)) stop("--contrast is required when --find_dmps true")
cts <- trimws(strsplit(contrast, ",", fixed = TRUE)[[1]])
cm <- limma::makeContrasts(contrasts = cts, levels = design)

# Fit the inferential model on M-values, and the identical design/contrast on
# beta values to obtain covariate-adjusted effect sizes on the interpretable scale.
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(m, design), cm))
beta_fit <- limma::contrasts.fit(limma::lmFit(b, design), cm)
fdr <- as.numeric(get("--fdr"))
dbthr <- as.numeric(get("--min-abs-delta-beta"))
write.table(data.frame(design = get("--design"), contrast = cts), file.path(out, "model_definitions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

group_summary <- function(result, contrast_name) {
  if (!"group" %in% names(ss)) return(result)
  terms <- strsplit(contrast_name, "-", fixed = TRUE)[[1]]
  group_columns <- paste0("group", make.names(as.character(unique(ss$group))))
  labels <- as.character(unique(ss$group))
  names(labels) <- group_columns
  if (length(terms) != 2 || !all(terms %in% names(labels))) return(result)
  for (term in terms) {
    label <- labels[[term]]
    ids <- ss$group == label
    suffix <- gsub("[^A-Za-z0-9._-]", "_", label)
    result[[paste0("mean_beta_", suffix)]] <- rowMeans(b[result$probe_id, ids, drop = FALSE], na.rm = TRUE)
    result[[paste0("n_", suffix)]] <- sum(ids)
  }
  result
}

for (ct in colnames(cm)) {
  z <- limma::topTable(fit, coef = ct, number = Inf, sort.by = "P")
  z$probe_id <- rownames(z)
  z$delta_beta <- beta_fit$coefficients[match(z$probe_id, rownames(beta_fit$coefficients)), ct]
  z <- group_summary(z, ct)
  safe <- gsub("[^A-Za-z0-9._-]", "_", ct)
  keep <- !is.na(z$adj.P.Val) & z$adj.P.Val <= fdr & !is.na(z$delta_beta) & abs(z$delta_beta) >= dbthr
  write.table(z, file.path(out, "complete", paste0(safe, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(z[keep, , drop = FALSE], file.path(out, "significant", paste0(safe, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
}
