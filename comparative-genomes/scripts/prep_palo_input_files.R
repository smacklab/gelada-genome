#!/usr/bin/env Rscript

i = commandArgs(trailingOnly=TRUE)[1]

to.do = scan(paste0('palo-tmp/',i,'/to_do.txt'),what='',sep='\n')
all.proteins = readRDS(paste0('palo-tmp/',i,'/proteins.rds'))

# n.total = nrow(read.table('assembly_list.txt',sep='\t',stringsAsFactors=FALSE))
orthogroups.totals = read.table('id_lists/orthogroups/sco_size.txt',sep='\t',header=FALSE,row.names=1)

n.total = orthogroups.totals[i,]

if (length(to.do) < n.total) {
	# PALO clips off ENSEMBL version numbers only, creating possible mismatches. So match proteins without version numbers
	palo.chosen = read.table(paste0('palo-tmp/',i,'/files/palo_combs.txt'),sep='\t',row.names=1,stringsAsFactors=FALSE)
	palo.chosen = gsub('\\.[0-9]*$','',as.character(unlist(palo.chosen[nrow(palo.chosen),])))
	# If number to do is less than the total, need to keep only the chosen proteins
	proteins.to.do = subset(all.proteins,gene.id %in% to.do)
	proteins.done = subset(all.proteins,gsub('\\.[0-9]*$','',protein.accession) %in% palo.chosen)
	proteins = rbind(proteins.done,proteins.to.do)
} else {
	# Otherwise
	proteins = all.proteins
}

# proteins now includes all proteins less the proteins that were excluded by previous runs of PALO

# Only write the first n taxa for which together there are fewer than 1e8 combinations
# Then write the remaining gene and protein data as the remaining chunk.
# For jobs that are incomplete, use a booster script to add more taxa until complete

# index is the number of genes to run
index = sum(unlist(lapply(seq(n.total),function(n) prod(table(proteins$gene.id)[1:n]))) < 1e8,na.rm=TRUE)

genes.to.do = levels(proteins$gene.id)[1:index]
genes.skipped = levels(proteins$gene.id)[-(1:index)]

x = droplevels(subset(proteins,gene.id %in% genes.to.do))

x$gene.id = as.character(x$gene.id)

write.table(
		x[c('gene.id','protein.accession','protein.length')],
		file = paste0('palo-tmp/',i,'/files/species.txt'),
		sep='\t',
		col.names=FALSE,
		row.names=FALSE,
		quote=FALSE)
y = unique(x[c('gene.id','orthogroup','assembly')])
y$assembly = as.character(y$assembly)
write.table(
		reshape2::acast(y,orthogroup~assembly,value.var='gene.id'),
		file = paste0('palo-tmp/',i,'/files/homologs.txt'),
		sep='\t',
		col.names=FALSE,
		row.names=FALSE,
		quote=FALSE)

write(genes.skipped,sep='\n',file=paste0('palo-tmp/',i,'/to_do_next.txt'))
