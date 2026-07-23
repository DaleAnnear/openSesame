#!/usr/bin/env Rscript

# Run with: Rscript tests/r/test_filter_probes.R bin/filter_probes.R
script <- commandArgs(trailingOnly = TRUE)[1L]
if (is.na(script) || !nzchar(script)) stop("Pass the path to filter_probes.R as the first argument.")

work <- tempfile("opensesame-filter-")
dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)

probe_ids <- paste0("cg", sprintf("%03d", seq_len(4L)))
sample_ids <- paste0("S", seq_len(4L))
beta <- matrix(0.5, nrow = 4L, ncol = 4L, dimnames = list(probe_ids, sample_ids))
detection <- matrix(c(
  0.01, 0.01, 0.01, 0.01,
  0.01, 0.10, 0.01, 0.10,
  0.01, 0.01, 0.01, 0.01,
  0.10, 0.10, 0.01, 0.01
), nrow = 4L, byrow = TRUE, dimnames = list(probe_ids, sample_ids))
beta["cg003", c("S1", "S2")] <- NA_real_
beta["cg004", c("S1", "S2")] <- NA_real_
analysis <- list(beta = beta, mvalue = beta, detection = detection)

analysis_path <- file.path(work, "analysis.rds")
output_path <- file.path(work, "out")
saveRDS(analysis, analysis_path)

result <- system2("Rscript", c(
  script,
  "--analysis", analysis_path,
  "--outdir", output_path,
  "--min-detection-rate", "0.75",
  "--max-missing-rate", "0.25",
  "--remove-sex", "false",
  "--remove-noncpg", "false"
), stdout = TRUE, stderr = TRUE)
if (!identical(attr(result, "status"), NULL)) stop(paste(result, collapse = "\n"))

removed <- read.delim(file.path(output_path, "probe_lists", "removed_probes.tsv"), check.names = FALSE)
expected <- data.frame(
  probe_id = c("cg002", "cg003", "cg004"),
  reason = c("detection_rate", "missingness", "detection_rate;missingness")
)
if (!isTRUE(all.equal(removed, expected, check.attributes = FALSE))) {
  stop("Probe filtering was not applied per probe with the expected reasons.")
}

filtered <- readRDS(file.path(output_path, "analysis_filtered.rds"))
if (!identical(rownames(filtered$beta), "cg001")) stop("Only the passing probe should remain after filtering.")
cat("Probe-level detection and missingness filtering test passed.\n")
