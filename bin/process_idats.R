#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(sesame))

# Define options
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="IDAT file prefix (e.g., path/to/sample1)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="beta_values.csv",
              help="Output CSV file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$prefix)){
  print_help(opt_parser)
  stop("Prefix argument must be supplied.\n", call.=FALSE)
}

# Run openSesame
# openSesame does the basic processing: detection p-value masking, background subtraction, dye bias correction
message("Processing IDAT prefix: ", opt$prefix)

# We use the standard openSesame pipeline.
betas <- openSesame(opt$prefix)

# Write to CSV
message("Writing output to: ", opt$output)
write.csv(betas, file=opt$output, row.names=TRUE)

message("Done.")
