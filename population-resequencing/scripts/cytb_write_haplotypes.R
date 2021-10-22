#!/usr/bin/env Rscript

cytb = read.delim('cytb_haplotypes.txt',header=FALSE)
cytb$V2 = toupper(cytb$V2)
cytb.split = split(cytb,cytb$V2)

haplogroup.membership = lapply(cytb.split,function(x) {
	out = unique(x$V1)[order(-grepl('^h',unique(x$V1)),unique(x$V1))]
	out
})

names(haplogroup.membership) = names(cytb.split) = haplotypes = as.character(unlist(lapply(haplogroup.membership,function(x) x[1])))

haplotypes = haplotypes[order(-grepl('^NC_',haplotypes),-grepl('^FIL',haplotypes),-grepl('^h',haplotypes),-(grepl('^SKR',haplotypes)|-grepl('^CHK',haplotypes)),-grepl('^GUA',haplotypes),haplotypes)]

cytb.split = cytb.split[haplotypes]
haplogroup.membership = haplogroup.membership[haplotypes]

haplotypes[!(grepl('^NC_',haplotypes) | grepl('^FIL',haplotypes) | grepl('^h',haplotypes))] = paste0('H',formatC(sum(grepl('^h',haplotypes)) + 1:sum(!(grepl('^NC_',haplotypes) | grepl('^FIL',haplotypes) | grepl('^h',haplotypes))),width=2,flag=0))

names(haplogroup.membership) = haplotypes

out = data.frame(haplotype=haplotypes,sequence=unlist(lapply(cytb.split,function(x) unique(x$V2))))
rownames(out) = NULL

write.table(out,file='cytb_haplotypes_only.txt',row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')