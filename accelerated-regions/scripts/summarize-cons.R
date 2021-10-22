#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

window.size = as.integer(arguments[1])
percent.intact = as.integer(arguments[2])
percent.identity = as.integer(arguments[3])

if (length(arguments) > 3) {
	n.cores = as.integer(arguments[4])
} else {
	n.cores = 24
}

source('config/settings.cfg')

blocks = list.files('results/cons/',pattern=paste0(maf_version,'\\..+?\\.window',window.size,'\\.intact',percent.intact,'\\.identity',percent.identity,'\\.[0-9]{6}\\.lrt.txt'))

library(parallel)

# Column names
#null_scale alt_scale alt_subscale lnlratio pval
results = do.call(rbind,mclapply(blocks,function(x) {
	this = read.table(file.path('results/cons',x),sep='\t',header=FALSE,stringsAsFactors=FALSE,col.names=c('block','chr','pos','str','block_size','lnl_null','lnl_alt'))
	this
},mc.cores=n.cores))

results$lnl_ratio = with(results,lnl_alt - lnl_null)
results$tval = with(results,2 * lnl_ratio)
results$pval = with(results,1 - pchisq(tval,1))

results$qval = with(results,p.adjust(pval,'fdr'))

results$coord = with(results,paste0(chr,':',pos,'-',pos+block_size-1))

results$maf = with(results,paste0(maf_version,'.',chr,'.window',window.size,'.intact',percent.intact,'.identity',percent.identity,'.',gsub('.+?_([0-9]{6})$','\\1',block),'.maf'))

dir.create('checkpoints',showWarnings=FALSE)

saveRDS(results,file=file.path('checkpoints',paste0('accelerated_regions_cons_',maf_version,'.window',window.size,'.intact',percent.intact,'.identity',percent.identity,'.rds')))
