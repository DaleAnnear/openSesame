#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE); get <- function(k) {i<-match(k,args);if(is.na(i)||i==length(args))stop('Missing ',k);args[i+1]}
if(!requireNamespace('sesame',quietly=TRUE)) stop('sesame is required; use the pinned container profile')
sample_id<-get('--sample-id'); red<-get('--red'); green<-get('--green'); prep<-get('--prep-code'); out<-get('--output'); metrics<-get('--metrics'); versions<-get('--versions')
allowed <- c('Q','C','D','P','B','H','T'); unknown <- setdiff(strsplit(prep,'')[[1]],allowed); if(length(unknown)) stop('Invalid sesame preparation code character(s): ',paste(unknown,collapse=','))
prefix <- sub('_(Red|red)\\.idat(\\.gz)?$','',basename(red)); file.copy(red,paste0(prefix,'_Red.idat',if(grepl('\\.gz$',red))'.gz','')); file.copy(green,paste0(prefix,'_Grn.idat',if(grepl('\\.gz$',green))'.gz',''))
sdf <- sesame::readIDATpair(prefix); sdf <- sesame::prepSesame(sdf, prep=prep); beta <- sesame::getBetas(sdf, mask=TRUE); mval <- log2(pmin(pmax(beta,1e-6),1-1e-6)/(1-pmin(pmax(beta,1e-6),1-1e-6)))
det <- if('pval' %in% colnames(sdf)) sdf$pval else rep(NA_real_,length(beta)); names(det)<-names(beta)
obj <- list(sample_id=sample_id,beta=beta,mvalue=mval,detection=det,prep_code=prep,platform_requested='EPICv2',sesame_version=as.character(utils::packageVersion('sesame')))
saveRDS(obj,out); qc<-data.frame(sample_id=sample_id,probe_count=length(beta),missing_beta_rate=mean(is.na(beta)),detection_call_rate=mean(is.na(det)|det<=0.05),prep_code=prep); write.table(qc,metrics,sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(package=c('sesame','R'),version=c(as.character(utils::packageVersion('sesame')),R.version.string)),versions,sep='\t',row.names=FALSE,quote=FALSE)
