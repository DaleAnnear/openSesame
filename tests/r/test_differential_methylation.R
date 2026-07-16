#!/usr/bin/env Rscript

# Run with: Rscript tests/r/test_differential_methylation.R bin/differential_methylation.R
script <- commandArgs(trailingOnly = TRUE)[1L]
if (is.na(script) || !nzchar(script)) stop("Pass the path to differential_methylation.R as the first argument.")
if (!requireNamespace("limma", quietly = TRUE)) stop("limma is required for this test.")

work <- tempfile("opensesame-dmp-")
dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)

set.seed(1)
probe_ids <- paste0("cg", sprintf("%03d", seq_len(20L)))
sample_ids <- paste0("S", seq_len(8L))
groups <- c(rep("Case", 4L), rep("Control", 4L))
batches <- rep(c("B1", "B2"), 4L)
beta <- matrix(runif(160L, 0.2, 0.7), nrow = 20L, dimnames = list(probe_ids, sample_ids))
beta[, groups == "Case"] <- pmin(0.99, beta[, groups == "Case"] + 0.12)
mvalue <- log2(beta / (1 - beta))
analysis <- list(mvalue = mvalue, beta = beta)

analysis_path <- file.path(work, "analysis.rds")
samplesheet_path <- file.path(work, "samplesheet.csv")
output_path <- file.path(work, "out")
saveRDS(analysis, analysis_path)
write.csv(data.frame(sample_id = sample_ids, group = groups, batch = batches), samplesheet_path, row.names = FALSE, quote = FALSE)

result <- system2("Rscript", c(
  script,
  "--implementation-version", "3",
  "--analysis", analysis_path,
  "--samplesheet", samplesheet_path,
  "--outdir", output_path,
  "--design", shQuote("~ 0 + group + batch"),
  "--contrast", "groupCase-groupControl",
  "--fdr", "0.05",
  "--min-abs-delta-beta", "0.05"
), stdout = TRUE, stderr = TRUE)
if (!identical(attr(result, "status"), NULL)) stop(paste(result, collapse = "\n"))

dmp <- read.delim(file.path(output_path, "complete", "groupCase-groupControl.tsv"), check.names = FALSE)
required <- c("delta_beta_adjusted", "delta_beta_observed", "mean_beta_Case", "mean_beta_Control")
if (!all(required %in% names(dmp))) stop("The expected DMP effect-size columns were not written.")
if (!all(is.finite(dmp$delta_beta_observed)) || any(abs(dmp$delta_beta_observed) > 1)) stop("Observed delta beta is not bounded to [-1, 1].")

expected <- rowMeans(beta[dmp$probe_id, groups == "Case", drop = FALSE]) - rowMeans(beta[dmp$probe_id, groups == "Control", drop = FALSE])
if (!isTRUE(all.equal(dmp$delta_beta_observed, unname(expected), tolerance = 1e-6))) stop("Observed delta beta does not equal the difference of the two group mean beta values.")
cat("DMP adjusted/observed delta-beta test passed.\n")
