#!/usr/bin/env Rscript

results.file = file.path('results','cafe_results_treefam.txt')

ignore.checkpoints = FALSE

target.genome = 'thegel'

# Read in CAFE results file
cafe = read.delim(results.file,skip=10,header=FALSE,stringsAsFactors=FALSE)

# Column headings:
# 'ID'	'Newick'	'Family-wide P-value'	'Viterbi P-values'	'cut P-value'	'Likelihood Ratio'

# Only need the first four columns
cafe = cafe[1:4]
names(cafe) = c('id','nwk','p.family','p.branch')

# Calculate node ID and find the position of its p-value
node.ids = gsub('^#.+?\\(','(',scan(results.file,quiet=TRUE,nlines=1,skip=2,what='',sep='\n'))
node.names = gsub(';$','',scan('data/tree_node_names.nwk',quiet=TRUE,sep='\n',what=''))
node.integers = gsub('[A-z]','',node.ids)

# Call useful functions

source('scripts/custom_functions.R')

target.node = find.target.node(target.genome)
parent.node = find.parent.node(target.node)
sister.node = find.sister.node(target.node)

library(parallel)

# Fetch the number of genes per family in gelada
cafe$target.count = unlist(mclapply(cafe$nwk,function(x) {
	get.node.count(x,target.node)
},mc.cores=16))
cafe$parent.count = unlist(mclapply(cafe$nwk,function(x) {
	get.node.count(x,parent.node)
},mc.cores=16))
cafe$sister.count = unlist(mclapply(cafe$nwk,function(x) {
	get.node.count(x,sister.node)
},mc.cores=16))

target.index = find.node.index(target.node)
parent.index = find.node.index(parent.node)
sister.index = find.node.index(sister.node)

# Fetch the p-value of the gelada branch (found with target.index)
cafe$target.p = unlist(mclapply(cafe$p.branch,function(x) {
	get.node.p(x,target.index)
},mc.cores=16))
cafe$parent.p = unlist(mclapply(cafe$p.branch,function(x) {
	get.node.p(x,parent.index)
},mc.cores=16))
cafe$sister.p = unlist(mclapply(cafe$p.branch,function(x) {
	get.node.p(x,sister.index)
},mc.cores=16))

cafe$target.change = factor(0,levels=c('expansion','contraction','no change'))

# If the count in target is higher than in the MRCA, classify as expansion
cafe$target.change[cafe$target.count > cafe$parent.count] = 'expansion'
# If the count in target is lower than in the MRCA, classify as contraction
cafe$target.change[cafe$target.count < cafe$parent.count] = 'contraction'
# If the count in target is equal to that in the MRCA, classify as no change
cafe$target.change[cafe$target.count == cafe$parent.count] = 'no change'

cafe$target.fdr = p.adjust(cafe$target.p,'fdr')
cafe$family.fdr = p.adjust(cafe$p.family,'fdr')

cafe$target.expansion.p = cafe$target.p
cafe$target.expansion.p[cafe$target.change == 'contraction'] = 1 - cafe$target.p[cafe$target.change == 'contraction']






all.node.ids = unlist(str_extract_all(node.ids,'[0-9]+'))
all.node.shortnames = gsub('<[0-9]+>','',unlist(str_extract_all(node.ids,'[a-z]*<[0-9]+>')))
all.node.shortnames[all.node.shortnames == ''] = 'ancestral'
all.node.names = unlist(str_extract_all(node.names,'[A-z/ -]+'))
all.node.indices = unlist(lapply(all.node.ids,find.node.index))
all.node.parents = unlist(lapply(all.node.ids,find.parent.node))

node.parents = data.frame(
		node.id = all.node.ids,
		node.name = all.node.shortnames,
		node.fullname = all.node.names,
		node.index = all.node.indices,
		parent.id = all.node.parents,
		stringsAsFactors=FALSE)

if (!file.exists('checkpoints/cafe_treefam_results_all_nodes.rds')) {
	cafe.lineage.results = cafe[c('id','p.family')]

	for (i in 1:nrow(node.parents)) {
		node.id = node.parents$node.id[i]
		parent.id = node.parents$parent.id[i]
		node.index = node.parents$node.index[i]
		if (node.parents$node.name[i] == 'ancestral') {
			node.name = paste0('ancestral',node.parents$node.id[i])
		} else {
			node.name = node.parents$node.name[i]
		}
		message('Now processing node ',node.name)

		if (!is.na(node.index)) {
			this.count = unlist(mclapply(cafe$nwk,function(x) {
				get.node.count(x,node.id)
			},mc.cores=16))
			parent.count = unlist(mclapply(cafe$nwk,function(x) {
				get.node.count(x,parent.id)
			},mc.cores=16))

			cafe.lineage.results[[paste0(node.name,'.change')]] = factor(0,levels=c('expansion','contraction','no change'))
			# If the count in target is higher than in the MRCA, classify as expansion
			cafe.lineage.results[[paste0(node.name,'.change')]][this.count > parent.count] = 'expansion'
			# If the count in target is lower than in the MRCA, classify as contraction
			cafe.lineage.results[[paste0(node.name,'.change')]][this.count < parent.count] = 'contraction'
			# If the count in target is equal to that in the MRCA, classify as no change
			cafe.lineage.results[[paste0(node.name,'.change')]][this.count == parent.count] = 'no change'
		
			cafe.lineage.results[[paste0(node.name,'.p')]] = unlist(mclapply(cafe$p.branch,function(x) {
				get.node.p(x,node.index)
			},mc.cores=16))
		}
	}
	saveRDS(cafe.lineage.results,file='checkpoints/cafe_treefam_results_all_nodes.rds')
} else {
	cafe.lineage.results = readRDS('checkpoints/cafe_treefam_results_all_nodes.rds')
}

library(biomaRt)

hsap = useEnsembl('ENSEMBL_MART_ENSEMBL','hsapiens_gene_ensembl',version=102)
gene2uniprot = getBM(attributes = c('ensembl_gene_id','uniprotswissprot'),mart=hsap)
gene2uniprot = subset(gene2uniprot,uniprotswissprot != '')


treefam2uniprot = read.delim('uniprotACC2treefam.txt',stringsAsFactors=FALSE)

treefam2uniprot.hsap = subset(treefam2uniprot,external_db_id %in% gene2uniprot$uniprotswissprot)

treefam2uniprot.hsap = merge(treefam2uniprot.hsap,gene2uniprot,by.x='external_db_id',by.y='uniprotswissprot')

# hsap2gelada = getBM(attributes=c('ensembl_gene_id','mmulatta_homolog_ensembl_gene'),filters=treefam2uniprot.hsap$ensembl_gene_id,mart=hsap)



tgel = useEnsembl('ENSEMBL_MART_ENSEMBL','tgelada_gene_ensembl',version=102)
treefam2tgel = read.delim('treefam_hmm_Tgel_1.0_bestfams.txt')
ens2refseq = getBM(attributes=c('ensembl_gene_id','refseq_peptide_predicted'),filters='refseq_peptide_predicted',values=gsub('\\.1','',treefam2tgel$seq_id),mart=tgel)
treefam2tgel$refseq_peptide_predicted = with(treefam2tgel,gsub('\\.[0-9]+$','',seq_id))


ens2treefam = subset(merge(ens2refseq,treefam2tgel,by='refseq_peptide_predicted'),select=c('hmm_name','ensembl_gene_id','seq_id','bit_score'))

cafe.ensembl = merge(cafe,ens2treefam,by.x='id',by.y='hmm_name',all.x=TRUE,all.y=FALSE)



hsap2go = getBM(attributes = c('ensembl_gene_id','go_id'),filters='ensembl_gene_id',values=unique(treefam2uniprot.hsap$ensembl_gene_id),mart=hsap)
hsap2go = subset(hsap2go,go_id != '')


treefam2go = merge(treefam2uniprot.hsap,hsap2go,by='ensembl_gene_id')

o2go = unique(treefam2go[c('gene_tree_stable_id','go_id')])
names(o2go)[[1]] = 'orthogroup'

# Leave this in as definition of expansion for the purposes of overrepresentation/enrichment analysis
cafe$expansion = cafe$target.fdr < 0.2 & cafe$target.change %in% 'expansion'

cafe.short = cafe[c('id','p.family','target.p','target.change')]
cafe.go = merge(cafe.short,o2go,by.x='id',by.y='orthogroup')

library(topGO)

ortho.go = cafe.go[,c('id','go_id')]
ortho.id.to.go = mclapply(unique(ortho.go$id),function(x) sort(ortho.go$go_id[ortho.go$id == x]),mc.cores=16)
names(ortho.id.to.go) = unique(ortho.go$id)

orthogroups = cafe$target.expansion.p
names(orthogroups) = cafe$id

# CHANGE CONTRACTIONS TO P = 1 for gene expansion analysis
go.data = new('topGOdata',description='Simple session',ontology='BP',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)


# Fetch names for GO terms identified as relations by topGO
if (ignore.checkpoints || !file.exists('checkpoints/treefam_go_names.rds')) {
	# Create dummy vector for topGO

	# Put together all relevant GO terms
	go.ids = sort(go.data@graph@nodes)

	# Exclude blank GO terms
	go.ids = go.ids[nchar(go.ids) == 10]

	# Grab names and namespaces from the GO database
	library(GO.db)

	go.terms = data.frame(
		go_id = go.ids,
		go_namespace = Ontology(go.ids),
		go_name = Term(go.ids),
		stringsAsFactors=FALSE
	)

	# Sometimes, there are names that are not in the GO database. Fetch these names from AmiGO

	library(XML)

	to.do = go.terms$go_id[!complete.cases(go.terms)]

	if (length(to.do)) {
		# Fetch missing GO metadata from AmiGO
		for (i in to.do) {
			message(i)
			search.url = paste0('http://amigo.geneontology.org/amigo/term/',i)

			amigo.data = xpathApply(htmlParse(search.url), '//dd', xmlValue)
			go.terms$go_namespace[match(i,go.terms$go_id)] = unlist(lapply(strsplit(amigo.data[[3]],'_'),function(x) paste(toupper(substr(x,1,1)),collapse='')))
			go.terms$go_name[match(i,go.terms$go_id)] = amigo.data[[2]]
		}
	}
	rownames(go.terms) = go.terms$go_id
	saveRDS(go.terms,file='checkpoints/treefam_go_names.rds')
} else {
	# If checkpoint is found, use it
	message('Checkpoint found!\nLoading GO metadata from file.')

	go.terms = readRDS('checkpoints/treefam_go_names.rds')
}


ksWeight01 = runTest(go.data,algorithm='weight01',statistic='ks')
# ksElim = runTest(go.data,algorithm='elim',statistic='ks')

orthogroups = as.integer(!cafe$expansion)
names(orthogroups) = cafe$id

# CHANGE CONTRACTIONS TO P = 1 for gene expansion analysis
go.data = new('topGOdata',description='Simple session',ontology='BP',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)

# FOCUS ON FET
# fetWeight01 = runTest(go.data,algorithm='weight01',statistic='fisher')
fetParentchild = runTest(go.data,algorithm='parentchild',statistic='fisher')


ortho.topgo = within(go.terms,{
#	fet.weight01 = fetWeight01@score[match(ortho.topgo$go_id,names(fetWeight01@score))]
	fet.parentchild = fetParentchild@score[match(go.terms$go_id,names(fetParentchild@score))]
	ks.weight01 = ksWeight01@score[match(go.terms$go_id,names(ksWeight01@score))]
#	ks.elim = ksElim@score[match(ortho.topgo$go_id,names(ksElim@score))]
})

ortho.topgo = within(ortho.topgo,{
#	fet.weight01.fdr = p.adjust(fet.weight01,'fdr')
	fet.parentchild.fdr = p.adjust(fet.parentchild,'fdr')
	ks.weight01.fdr = p.adjust(ks.weight01,'fdr')
#	ks.elim.fdr = p.adjust(ks.elim,'fdr')
})


ortho.topgo = ortho.topgo[order(ortho.topgo$ks.weight01.fdr),]


cafe = cafe[order(cafe$target.expansion.p),]


treefam2gene = treefam2uniprot.hsap

treefam2gene = getBM(attributes = c('ensembl_gene_id','external_gene_name'),filters='ensembl_gene_id',values=unique(treefam2uniprot.hsap$ensembl_gene_id),mart=hsap)

treefam2gene = merge(treefam2uniprot.hsap,treefam2gene,by='ensembl_gene_id')

treefam2.gene.join = unlist(lapply(split(treefam2gene,treefam2gene$gene_tree_stable_id),function(x) sort(paste(unique(x$external_gene_name),collapse=', '))))


cafe$external_gene_name = treefam2.gene.join[cafe$id]

saveRDS(ortho.topgo,file='checkpoints/cafe_go_results.rds')

saveRDS(cafe,file='checkpoints/cafe_results.rds')
