#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

genome = arguments[1]
fai = arguments[2]

chromosomes = read.delim(fai,header=FALSE,stringsAsFactors=FALSE,col.names=c('chromosome','size',paste0('V',3:5)))[1:2]

chromosomes$int = formatC(1:nrow(chromosomes),width=2,flag='0')

# Exclude mitochondria and small scaffolds
chromosomes = subset(chromosomes,size > 1e7)

regions = unlist(lapply(1:nrow(chromosomes),function(i) {
	this.chr = chromosomes$chromosome[i]
	this.size = chromosomes$size[i]
	this.start = seq(1,this.size,3e6)
	this.end = seq(3e6,this.size + 3e6,3e6)
	this.end[this.end > this.size] = this.size
	paste0(this.chr,':',gsub(' ','',format(this.start,scientific=FALSE)),'-',gsub(' ','',format(this.end,scientific=FALSE)))
}))

write(regions,sep='\n',file=paste0('data/',genome,'_regions.txt'))
