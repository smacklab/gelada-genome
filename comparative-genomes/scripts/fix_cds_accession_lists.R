#!/usr/bin/env Rscript

# Read in single-copy orthologs (SCOs)
# sco = scan(file='id_lists/orthogroups/sco_all.txt',sep='\n',what='')
sco = scan(file='id_lists/orthogroups/sco_most.txt',sep='\n',what='')

library(parallel)

system('mkdir -p id_lists/cds_palo_isoform_corrected')

dev.null = mclapply(sco,function(x) {
	cds.accessions = scan(paste0('id_lists/cds_palo_isoform/',x,'.txt'),what='\n',quiet=TRUE)
	# Refseq
	r = which(!grepl('^ENS',cds.accessions))
	# Ensembl
	e = grep('^ENS',cds.accessions)
	cds.corrected = system(paste0('grep ">" tmp/cds/',x,'.fa | cut -d " " -f 1 | sed "s/>//g"'),intern=TRUE)
	out = character(length(cds.accessions))
	out[r] = cds.corrected
	out[e] = cds.accessions[e]
	write(out,file=paste0('id_lists/cds_palo_isoform_corrected/',x,'.txt'),sep='\n')
	length(out[as.logical(nchar(out))])
},mc.cores=16)
