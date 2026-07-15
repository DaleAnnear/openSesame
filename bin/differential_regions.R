#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);get<-function(k){i<-match(k,args);if(is.na(i))stop('Missing ',k);args[i+1]};out<-get('--outdir');dir.create(out,recursive=TRUE,showWarnings=FALSE)
need<-c('DMRcate','IlluminaHumanMethylationEPICv2anno.20a1.hg38','EPICv2manifest');missing<-need[!vapply(need,requireNamespace,logical(1),quietly=TRUE)];if(length(missing))stop('EPICv2 DMR analysis requires pinned packages missing from this image: ',paste(missing,collapse=', '),'. No EPICv1 coordinate fallback is permitted.')
stop('DMRcate EPICv2 execution is gated pending a tested annotation-to-DMP adapter; see docs/feature_parity.md. This explicit stop prevents unvalidated region calls.')
