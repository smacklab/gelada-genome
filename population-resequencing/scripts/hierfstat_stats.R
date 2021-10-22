#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

prefix = arguments[1]
this.i = as.integer(arguments[2])
this.step = arguments[3]

vcf.file = file.path('vcf_snps',list.files('vcf_snps',pattern=paste0(prefix,'.*.good.step',this.step,'.vcf.gz$')))[this.i]

if (!file.exists(file.path(
	'fst_results',
	gsub('.vcf.gz$','.fst.rds',gsub('.+?/','',vcf.file))
))) {

# Read in options to ensure that unwanted chromosomes are excluded
source('scripts/_include_options.R')

# vcf.file = 'vcf_final/gelada.panubis1.bootstrap.all.step2.vcf.gz'
# vcf.file = 'vcf/gelada.tgel1.bootstrap.region.0001.chr01.good.step2.vcf.gz'

library(gaston)

vcf = read.vcf(vcf.file,convert.chr=FALSE)

# library(adegenet)

# Gelada is 1, baboon is 2
# ifelse(substr(vcf@ped$id,1,3) %in% c('FIL'),2,1)

snp.data = data.frame(pop=ifelse(substr(vcf@ped$id,1,3) %in% c('FIL'),2,1),as.data.frame(as.matrix(vcf)[,!vcf@snps$chr %in% exclude.chrs]))

library(hierfstat)

# vcf.stats = basic.stats(snp.data)

wc.fst = pairwise.WCfst(snp.data)

chr = as.integer(gsub('vcf_snps/.*\\.chr([0-9]{2})\\..*','\\1',vcf.file))
n.snps = ncol(snp.data) - 1
fst = as.numeric(na.omit(unique(as.numeric(wc.fst))))

out = data.frame(vcf=vcf.file,chr,fst,n.snps)

dir.create('fst_results',showWarnings=FALSE)

saveRDS(out,file=file.path(
	'fst_results',
	gsub('.vcf.gz$','.fst.rds',gsub('.+?/','',vcf.file))
))

}

#           1         2
# 1        NA 0.3240103
# 2 0.3240103        NA