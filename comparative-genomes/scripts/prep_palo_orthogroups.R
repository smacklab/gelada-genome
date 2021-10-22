#!/usr/bin/env Rscript

# Read in orthogroups (orthogroup IDs and associated proteins)
orthogroups = read.delim('id_lists/orthogroups/orthogroups.tsv',stringsAsFactors=FALSE)

# Read in single-copy orthologs (SCOs)
sco = scan(file='id_lists/orthogroups/sco.txt',sep='\n',what='')

# Now calculate SCO+ (single-copy orthologs with missingness in genomes)
orthogroups.matrix = do.call(cbind,lapply(orthogroups[2:ncol(orthogroups)],function(o) {
	o[!nchar(o)] = NA
	unlist(lapply(strsplit(o,','),function(x) length(x[!is.na(x)])))
}))
rownames(orthogroups.matrix) = orthogroups[[1]]

sco.all = orthogroups.matrix[apply(orthogroups.matrix,1,function(x) all(x <= 1)),]
sco.plus = sco.all[rowSums(sco.all) >= ncol(sco.all) / 2 & rowSums(sco.all) != ncol(sco.all),]
sco.plus = sco.plus[order(rowSums(sco.plus),decreasing=TRUE),]

# Now calculate SCO? (single-copy orthologs but for duplications in a few)
sco.most = orthogroups.matrix[apply(orthogroups.matrix,1,function(x) sum(x == 1) >= ncol(orthogroups.matrix)/2) & !apply(orthogroups.matrix,1,function(x) all(x <= 1)),]

# Force the number to be 0
sco.most[sco.most > 1] = 0

# No point including the orthogroup if geladas don't count
sco.most = sco.most[as.logical(sco.most[,'GCF_003255815.1_Tgel_1.0_protein']),]

# Update all SCOs to include sco.plus and sco.most
sco.all = rbind(
	sco.all[sco,],
	sco.all[rownames(sco.plus),],
	sco.most
)

write(rownames(sco.plus),file='id_lists/orthogroups/sco_plus.txt',sep='\n')
write(rownames(sco.most),file='id_lists/orthogroups/sco_most.txt',sep='\n')
write(rownames(sco.all),file='id_lists/orthogroups/sco_all.txt',sep='\n')

write.table(reshape2::melt(rowSums(sco.all)),file='id_lists/orthogroups/sco_size.txt',sep='\t',col.names=FALSE,quote=FALSE,row.names=TRUE)

# For now, merge SCOs with missingness into sco vector
sco = c(sco,rownames(sco.plus),rownames(sco.most))

# Drop orthogroups that are not SCOs
row.names(orthogroups) = orthogroups[[1]]
sco.orthogroups = orthogroups[sco,]

# Erase proteins that have commas (indicating multiple copies)
sco.orthogroups = do.call(data.frame,lapply(sco.orthogroups,function(x) gsub('^[A-Z0-9_.]+,.*','',x)))
rownames(sco.orthogroups) = sco

sco.orthogroups.melted = reshape2::melt(sco.orthogroups,id='Orthogroup')
names(sco.orthogroups.melted) = c('orthogroup','assembly','protein.accession')

# There is a bug in which special characters (e.g., hyphens) get substituted
#     when converted to R data frame headings. As this is less likely to
#     happen to the assembly ID, match based on the genome accession
#     (both NCBI and ensembl should have at least one underscore in the real ID
#     so grabbing all characters before the SECOND underscore should be sufficient)

# Read in proper assembly metadata
assembly.list = read.table('assembly_list.txt',sep='\t',stringsAsFactors=FALSE)
assembly.list$assembly = gsub(' ','_',with(assembly.list,paste(V1,V2,sep='_')))
rownames(assembly.list) = gsub('(^.+?_.+?)_.+','\\1',gsub(' ','_',assembly.list$V1))

# Format assembly data in melted data frame to match
sco.orthogroups.melted$assembly = gsub('_protein$','',sco.orthogroups.melted$assembly)
sco.orthogroups.melted$assembly = assembly.list[gsub('(^.+?_.+?)_.+','\\1',sco.orthogroups.melted$assembly),]$assembly

sco.orthogroups.split = split(sco.orthogroups.melted,sco.orthogroups.melted$assembly)

library(parallel)

all.o2p = do.call(rbind,mclapply(names(sco.orthogroups.split),function(i) {
	x = sco.orthogroups.split[[i]]
	p2g = read.table(paste0('id_lists/feature_joins/',i,'_prot2gene.txt'),sep='\t',stringsAsFactors=FALSE,col.names=c('gene.id','protein.accession'))
	p2c = read.table(paste0('id_lists/feature_joins/',i,'_prot2cds.txt'),sep='\t',stringsAsFactors=FALSE,col.names=c('protein.accession','cds.id'))
	o2g = merge(x,p2g)
	o2g$protein.accession = NULL
	o2p = merge(o2g,p2g,by='gene.id')
	o2p = merge(o2p,p2c,by='protein.accession')
	fai = read.table(paste0('fasta/aa/',i,'_protein.faa.fai'),sep='\t',stringsAsFactors=FALSE,col.names=c('protein.accession','protein.length','start','c1','c2'))[1:2]
	o2p = merge(o2p,fai,by='protein.accession')
	o2p$gene.id = as.character(o2p$gene.id)
	o2p
},mc.cores=16))
o2p.split.gene = split(all.o2p,all.o2p$gene.id)

# For each gene, if there are multiple isoforms of the same length, just choose the first one to lighten the load
o2p.dedup = do.call(rbind,mclapply(o2p.split.gene,function(x) {
	x[as.integer(tapply(1:nrow(x),x$protein.length,function(x) x[1])),]
},mc.cores=4))

o2p.split.orthogroup = split(o2p.dedup,o2p.dedup$orthogroup)

# Palo requires 2 files

# species.txt is a tab-separated table with one row per isoform and three columns: the gene identifier, protein identifier, and the length
# homologs.txt has one row per orthogroup. For SCOs, there should be one column per species. Columns will be populated by Gene IDs from that species

dev.null = mclapply(names(o2p.split.orthogroup),function(i) {
	x = o2p.split.orthogroup[[i]]
	x$gene.id = factor(x$gene.id,levels=names(sort(table(x$gene.id))))

	system(paste0('mkdir -p palo-tmp/',i))
	system(paste0('mkdir -p palo-tmp/',i,'/files'))

	saveRDS(x,file=paste0('palo-tmp/',i,'/proteins.rds'))
	write(levels(x$gene.id),sep='\n',file=paste0('palo-tmp/',i,'/to_do.txt'))
	nrow(x)
},mc.cores=16)
