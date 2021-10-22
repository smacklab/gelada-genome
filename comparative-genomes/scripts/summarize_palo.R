#!/usr/bin/env Rscript

i = commandArgs(trailingOnly=TRUE)[1]

all.proteins = readRDS(paste0('palo-tmp/',i,'/proteins.rds'))
palo.chosen = read.table(paste0('palo-tmp/',i,'/files/palo_combs.txt'),sep='\t',row.names=1,stringsAsFactors=FALSE)
palo.chosen = gsub('\\.[0-9]*$','',as.character(unlist(palo.chosen[nrow(palo.chosen),])))

proteins = subset(all.proteins,gsub('\\.[0-9]*$','',protein.accession) %in% palo.chosen)

proteins = proteins[order(proteins$assembly),]

system('mkdir -p id_lists/aa_palo_isoform')
system('mkdir -p id_lists/cds_palo_isoform')

write(proteins$protein.accession,file=paste0('id_lists/aa_palo_isoform/',i,'.txt'),sep='\n')
write(proteins$cds.id,file=paste0('id_lists/cds_palo_isoform/',i,'.txt'),sep='\n')
