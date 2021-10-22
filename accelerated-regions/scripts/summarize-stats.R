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

blocks = list.files('reports/fin',pattern=paste0(maf_version,'.*.window',window.size,'.intact',percent.intact,'.identity',percent.identity,'.*.stats.txt'))

library(parallel)

# Column names
header.names = scan(file.path('reports/fin',blocks[1]),nlines=1,what='')

#null_scale alt_scale alt_subscale lnlratio pval
stats = do.call(rbind,mclapply(blocks,function(x) {
	this = read.table(file.path('reports/fin',x),sep='\t',header=FALSE,stringsAsFactors=FALSE,skip=1,col.names=header.names)
	this$block = gsub('\\.stats\\.txt$','',x)
	this$maf = gsub('\\.stats\\.txt$','.maf',x)
	this
},mc.cores=n.cores))

ref.blocks = list.files('reports/ref',pattern=paste0(maf_version,'.*.window',window.size,'.intact',percent.intact,'.identity',percent.identity,'.*.stats.',ref_species,'.txt'))

human = do.call(rbind,mclapply(ref.blocks,function(x) {
	this = read.table(file.path('reports/ref',x),sep='\t',header=FALSE,stringsAsFactors=FALSE,skip=1,col.names=header.names)
	this$block = gsub(paste0('\\.stats\\.',ref_species,'\\.txt$'),'',x)
	this$maf = gsub(paste0('\\.stats\\.',ref_species,'\\.txt$'),'.maf',x)
	this
},mc.cores=n.cores))

rownames(human) = human$maf
human = human[stats$maf,]
rownames(human) = NULL
human = subset(human,select=c('Chr','Start','Stop'))
names(human) = paste0(ref_species,'.',names(human))

stats = data.frame(stats,human)

saveRDS(stats,file=file.path('checkpoints',paste0('accelerated_regions_stats_',maf_version,'.window',window.size,'.intact',percent.intact,'.identity',percent.identity,'.rds')))
