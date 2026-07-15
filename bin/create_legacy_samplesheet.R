#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE); i <- match('--output',args); if(is.na(i)) stop('Use --output'); out <- args[i+1]; files <- args[-c(i,i+1)]
files <- normalizePath(files,mustWork=TRUE); red <- files[grepl('_(Red|red)\\.idat(\\.gz)?$',files)]; green <- files[grepl('_(Grn|Green|grn|green)\\.idat(\\.gz)?$',files)]
prefix <- function(z) sub('_(Red|red|Grn|Green|grn|green)\\.idat(\\.gz)?$','',z); keys <- sort(unique(c(prefix(red),prefix(green))))
x <- data.frame(sample_id=basename(keys), idat_red=red[match(keys,prefix(red))], idat_green=green[match(keys,prefix(green))], stringsAsFactors=FALSE)
if(anyNA(x$idat_red)|anyNA(x$idat_green)) stop('Glob contains unpaired IDAT files'); write.csv(x,out,row.names=FALSE)
