#!/usr/bin/env Rscript

library(parallel)

assemblies = read.table('assembly_list.txt',sep='\t',stringsAsFactors=FALSE,na.strings=c(''))
orthogroups = read.delim('id_lists/orthogroups/orthogroups.tsv',stringsAsFactors=FALSE,na.strings='')

assembly.names = assemblies$V4
names(assembly.names) = with(assemblies,gsub('-','.',gsub(' ','_',paste0(V1,'_',V2,'_protein'))))

names(orthogroups)[names(orthogroups) %in% names(assembly.names)] = assembly.names[names(orthogroups)[names(orthogroups) %in% names(assembly.names)]]

# Orthogroup to protein lists
o2p = mclapply(assembly.names,function(x) {
	out = lapply(orthogroups[[x]],function(y) {
		if (sum(!is.na(y))) {
			unlist(strsplit(gsub(' ','',y),','))
		} else {
			NA
		}
	})
	names(out) = orthogroups$Orthogroup
	out
},mc.cores=16)
names(o2p) = assembly.names
# Taxon to gene lists
t2g = mclapply(assembly.names,function(x) {
	assembly.data = subset(assemblies,V4 == x)
	p2g = read.table(with(assembly.data,gsub(' ','_',paste0('id_lists/feature_joins/',V1,'_',V2,'_prot2gene.txt'))),stringsAsFactors=FALSE)
	gene.lookup = as.character(p2g$V1)
	names(gene.lookup) = p2g$V2
	gene.lookup
},mc.cores=16)
names(t2g) = assembly.names

# Now make orthogroup to gene lists
o2g = mclapply(assembly.names,function(x) {
	taxon.genes = t2g[[x]]
	out = lapply(o2p[[x]],function(y) {
		if (sum(!is.na(y))) {
			unique(as.character(taxon.genes[y]))
		} else {
			NA
		}
	})
	names(out) = orthogroups$Orthogroup
	out
},mc.cores=16)
names(o2g) = assembly.names

saveRDS(o2g,file='checkpoints/orthogroups_to_genes.rds')

gene.lists = mclapply(assembly.names,function(x) {
	out = unlist(o2g[[x]])
	unique(as.character(out[!is.na(out)]))
},mc.cores=16)
names(gene.lists) = assembly.names

library(biomaRt)

if (!file.exists('checkpoints/biomart_go_all.rds')) {
	for (i in assembly.names[!is.na(assemblies$V7)]) {
		cat(i,'\n')
		assembly.info = subset(assemblies,V4 == i)
		genes = gene.lists[[i]]
		mart = useEnsembl(biomart='ENSEMBL_MART_ENSEMBL',dataset=assembly.info$V7,version=102)
		if (assembly.info$V3 == 'refseq') {
			go = getBM(attributes=c('ensembl_gene_id','external_gene_name','entrezgene_id','go_id','name_1006','namespace_1003'),filters='entrezgene_id',values=genes,mart=mart)
			go$gene_id = go$entrezgene_id
		} else {
			go = getBM(attributes=c('ensembl_gene_id','external_gene_name','go_id','name_1006','namespace_1003'),filters='ensembl_gene_id',values=genes,mart=mart)
			go$entrezgene_id = NA
			go$gene_id = go$ensembl_gene_id
		}
		cat(nrow(go),'rows\n')
		if (nrow(go)) {
			go$taxon = i
			go = go[c('ensembl_gene_id','external_gene_name','entrezgene_id','go_id','name_1006','namespace_1003','taxon','gene_id')]
			saveRDS(go,file=paste0('checkpoints/',i,'_biomart_go.rds'))
		}
	}

	all.go = do.call(rbind,mclapply(assembly.names,function(x) {
		if (file.exists(paste0('checkpoints/',x,'_biomart_go.rds'))) {
			readRDS(paste0('checkpoints/',x,'_biomart_go.rds'))
		} else {
			data.frame(
				ensembl_gene_id = character(0),
				external_gene_name = character(0),
				entrezgene_id = logical(0),
				go_id = character(0),
				name_1006 = character(0),
				taxon = character(0))
		}
	},mc.cores=16))
	rownames(all.go) = NULL

	saveRDS(all.go,file='checkpoints/biomart_go_all.rds')
} else {
	all.go = readRDS('checkpoints/biomart_go_all.rds')
}
library(reshape2)

o2gene = do.call(rbind,mclapply(assembly.names,function(x) {
	out = melt(o2g[[x]])
	names(out) = c('gene_id','orthogroup')
	out$taxon = x
	out = out[c('orthogroup','gene_id','taxon')]
	out
},mc.cores=16))
rownames(o2gene) = NULL

o2go.all = merge(o2gene,subset(all.go,select=c('taxon','gene_id','go_id','name_1006','namespace_1003')),by=c('taxon','gene_id'))

o2go.all = subset(o2go.all,!go_id %in% c('',NA))

o2go = unique(subset(o2go.all,select=c('orthogroup','go_id','name_1006','namespace_1003')))

saveRDS(o2go,file='checkpoints/orthogroups_to_go.rds')

# Count number of proteins for each assembly
o2p.counts = do.call(cbind,lapply(o2p,function(x) {
	unlist(lapply(x,function(y) {
		length(y[!is.na(y)])
	}))
}))

# Filter to only SCOs (with missingness allowed)
o2p.sco = o2p.counts[apply(o2p.counts,1,function(x) all(x <= 1)),]

source('scripts/custom_functions.R')

# orthogroups with all but one taxon
o2p.sco.abm = o2p.sco[rowSums(o2p.sco) == ncol(o2p.sco) - 1,]

o2p.sco.abm.df = data.frame(
	orthogroup = rownames(o2p.sco.abm),
	missing.taxon = colnames(o2p.sco.abm)[apply(o2p.sco.abm,1,which.min)],
	stringsAsFactors=FALSE
)

write.table(o2p.sco.abm.df,file='id_lists/orthogroups/sco_abm.txt',sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)

nwk = scan('data/tree.nwk',what='')

for (i in assemblies$V4) {
	nwk.out = prune.branch(i,nwk)
	nwk.paml = gsub('(thegel:[0-9.]+)','\\1#1',nwk.out)
	write(nwk.paml,file=paste0('data/tree_abm_',i,'_paml.txt'))
	nwk.busted = gsub('(thegel)(:[0-9.]+)','\\1{Foreground}\\2',nwk.out)
	write(nwk.busted,file=paste0('data/tree_abm_',i,'_busted.txt'))
}
