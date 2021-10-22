e#!/usr/bin/env Rscript

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

blocks = list.files('results/phylop',pattern=paste0(maf_version,'\\..+?\\.window',window.size,'\\.intact',percent.intact,'\\.identity',percent.identity,'\\.[0-9]{6}\\.wig'))

library(parallel)

# Column names
#null_scale alt_scale alt_subscale lnlratio pval
results = do.call(rbind,mclapply(blocks,function(x) {
	chr = gsub(paste0(maf_version,'\\.(.+?)\\.window',window.size,'\\.intact',percent.intact,'\\.identity',percent.identity,'\\.([0-9]{6})\\.wig'),'\\1',x)
	int = gsub(paste0(maf_version,'\\.(.+?)\\.window',window.size,'\\.intact',percent.intact,'\\.identity',percent.identity,'\\.([0-9]{6})\\.wig'),'\\2',x)
	path = file.path('results/phylop',x)
	header = scan(path,what='',sep='\n',nmax=2,quiet=TRUE)[2]
	if (is.na(header)) {
		out = data.frame(block = character(0), chr = character(0), pos = integer(0), null_scale=numeric(0), alt_scale=numeric(0),alt_subscale=numeric(0),lnlratio=numeric(0),pval=numeric(0))
	} else {
		pos = as.integer(gsub('fixedStep chrom=.+?start=([0-9]+) step=1','\\1',header))
		out = read.table(path,skip=2,sep='\t',col.names=c('null_scale','alt_scale','alt_subscale','lnlratio','pval'))
		out = data.frame(block = paste(chr,int,sep='_'), chr = chr, pos = pos:(pos + nrow(out) - 1), out)
	}
	out
},mc.cores=n.cores))

dir.create('checkpoints',showWarnings=FALSE)

saveRDS(results,file=file.path('checkpoints',paste0('accelerated_regions_phylop_',maf_version,'.window',window.size,'.intact',percent.intact,'.identity',percent.identity,'.rds')))
