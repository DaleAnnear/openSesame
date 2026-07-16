#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get <- function(k) { i <- match(k, args); if (is.na(i) || i == length(args)) stop("Missing ", k); args[i + 1] }
if (!requireNamespace("limma", quietly = TRUE)) stop("limma is required for DMP analysis")

a <- readRDS(get("--analysis")); ss <- read.csv(get("--samplesheet"), check.names = FALSE, stringsAsFactors = FALSE); out <- get("--outdir")
dir.create(file.path(out, "complete"), recursive = TRUE, showWarnings = FALSE); dir.create(file.path(out, "significant"), recursive = TRUE, showWarnings = FALSE)
assay <- function(name) if (inherits(a, "SummarizedExperiment")) SummarizedExperiment::assay(a, name) else a[[name]]
m <- assay("mvalue"); b <- assay("beta")
if (!nrow(m)) stop("No probes remain after filtering. Inspect filtering/filtering_summary.tsv and filtering/probe_lists/removed_probes.tsv; for small cohorts, consider scientifically justified relaxation of --min_detection_rate and/or --max_missing_rate.")
if (ncol(m) < 2) stop("DMP analysis requires at least two samples after QC.")
if (!identical(rownames(m), rownames(b)) || !identical(colnames(m), colnames(b))) stop("M-value and beta matrices are not aligned.")
ss <- ss[match(colnames(m), ss$sample_id), , drop = FALSE]
if (anyNA(ss$sample_id)) stop("Samplesheet does not match the analysis matrix.")
design <- model.matrix(as.formula(get("--design")), ss)
if (qr(design)$rank < ncol(design)) stop("Design matrix is rank-deficient; revise design/covariates")
contrast <- get("--contrast"); if (contrast == "" || is.na(contrast)) stop("--contrast is required when --find_dmps true")
cts <- trimws(strsplit(contrast, ",", fixed = TRUE)[[1]]); cm <- limma::makeContrasts(contrasts = cts, levels = design)

# Statistical inference is on M-values. The paired beta-scale fit uses the same
# design and contrast to report a covariate-adjusted estimate. Ordinary least
# squares is unconstrained, so that estimate can lie outside [-1, 1] and must not
# be interpreted as a literal difference between observed beta values.
fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(m, design), cm)); beta_fit <- limma::contrasts.fit(limma::lmFit(b, design), cm)
fdr <- as.numeric(get("--fdr")); dbthr <- as.numeric(get("--min-abs-delta-beta"))
write.table(data.frame(design = get("--design"), contrast = cts), file.path(out, "model_definitions.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

group_contrast <- function(contrast_name) {
  if (!"group" %in% names(ss)) return(NULL)
  terms <- strsplit(contrast_name, "-", fixed = TRUE)[[1]]
  labels <- as.character(unique(ss$group))
  names(labels) <- paste0("group", make.names(labels))
  if (length(terms) != 2L || !all(terms %in% names(labels))) return(NULL)
  unname(labels[terms])
}

add_observed_group_delta <- function(result, contrast_name) {
  labels <- group_contrast(contrast_name)
  result$delta_beta_observed <- NA_real_
  if (is.null(labels)) return(result)

  means <- lapply(labels, function(label) {
    ids <- ss$group == label
    rowMeans(b[result$probe_id, ids, drop = FALSE], na.rm = TRUE)
  })
  names(means) <- labels
  for (label in labels) {
    suffix <- gsub("[^A-Za-z0-9._-]", "_", label)
    result[[paste0("mean_beta_", suffix)]] <- means[[label]]
    result[[paste0("n_", suffix)]] <- sum(ss$group == label)
  }
  # The contrast is numerator minus denominator, matching limma's contrast
  # notation. This is a difference of observed beta means and is bounded to
  # [-1, 1] whenever both group means are available.
  result$delta_beta_observed <- means[[1L]] - means[[2L]]
  result
}

summary_rows <- list()
for (ct in colnames(cm)) {
  # limma permits missing observations, but a row without an estimable test must
  # never be represented as a fabricated all-NA DMP result.
  testable <- is.finite(fit$t[, ct]) & is.finite(fit$p.value[, ct]) & is.finite(beta_fit$coefficients[, ct])
  summary_rows[[ct]] <- data.frame(contrast = ct, input_probes = nrow(m), testable_probes = sum(testable), excluded_nonestimable = sum(!testable), stringsAsFactors = FALSE)
  if (!any(testable)) stop("No probes have estimable statistics for contrast ", ct, ". This usually indicates excessive masking or insufficient non-missing values in one or more design groups; inspect detection metrics and relax filters only with scientific justification.")
  z <- limma::topTable(fit[testable, ], coef = ct, number = Inf, sort.by = "P")
  z$probe_id <- rownames(z)
  z$delta_beta_adjusted <- beta_fit$coefficients[match(z$probe_id, rownames(beta_fit$coefficients)), ct]
  z <- add_observed_group_delta(z, ct)
  safe <- gsub("[^A-Za-z0-9._-]", "_", ct)
  if (all(is.na(z$delta_beta_observed))) {
    stop("Cannot calculate a bounded observed delta beta for contrast ", ct,
         ". Use a two-group contrast of the form group<case>-group<reference> ",
         "with group labels matching the samplesheet.")
  }
  keep <- !is.na(z$adj.P.Val) & z$adj.P.Val <= fdr & !is.na(z$delta_beta_observed) & abs(z$delta_beta_observed) >= dbthr
  write.table(z, file.path(out, "complete", paste0(safe, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(z[keep, , drop = FALSE], file.path(out, "significant", paste0(safe, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
}
write.table(do.call(rbind, summary_rows), file.path(out, "dmp_testability_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
