#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

paml.results.file = file.path('results','paml_results_all.txt')
busted.results.file = file.path('results','busted_results_all.txt')
relax.results.file = file.path('results','relax_results_all.txt')

# Read in PAML results file
paml = read.delim(paml.results.file,header=FALSE,stringsAsFactors=FALSE)
busted = read.delim(busted.results.file,header=FALSE,stringsAsFactors=FALSE)
relax = read.delim(relax.results.file,header=FALSE,stringsAsFactors=FALSE)

names(paml) = c('id','h0','h1','lnl','p')
names(busted) = c('id','lnl','p')
names(relax) = c('id','p')

paml$p[!is.na(paml$lnl) & paml$lnl < 0] = 1

# Redo p

paml$p = with(paml,pchisq(2*(h1-h0),df=1,lower.tail=FALSE))

# Backup with all orthogroups
paml.all = paml

# Drop orthogroups with NA results
paml = subset(paml,!is.na(paml$p))

all.selection = merge(merge(paml,busted,by='id'),relax,by='id')
names(all.selection) = c('id','h0','h1','lnl','paml.p','busted.lnl','busted.p','relax.p')

# Call useful functions

source('scripts/custom_functions.R')

orthogroup.dims = read.table('results/orthogroups_dims.txt',sep='\t',fill=TRUE)
names(orthogroup.dims) = c('id','nchar','ntax')

if (!file.exists('sco_stats.txt')) {
	# More orthogroup info
	all.ortho = read.table('id_lists/orthogroups/orthogroups.tsv',sep='\t',stringsAsFactors=FALSE)
	rownames(all.ortho) = all.ortho$V1
	all.ortho$V1 = NULL
	names(all.ortho) = system('cut -f 4 assembly_list.txt',intern=TRUE)

	sco.ortho = all.ortho[scan('id_lists/orthogroups/sco_all.txt',what=''),]
	sco.ortho.copies = apply(sco.ortho,2,function(x) ifelse(stringr::str_length(x),stringr::str_count(x,pattern=',') + 1,0))
	sco.stats = data.frame(orthogroup = rownames(sco.ortho))

	sco.stats = within(sco.stats,{
		strict.sco = apply(sco.ortho.copies,1,function(x) all(x <= 1))
		n.copies = rowSums(sco.ortho.copies)
		n.missing = apply(sco.ortho.copies,1,function(x) sum(x == 0))
		n.duplicated = apply(sco.ortho.copies,1,function(x) sum(x > 1))
		max.duplicated = apply(sco.ortho.copies,1,function(x) if (any(x>1)) max(x[x>1]) else 1)
	})
	rownames(sco.stats) = sco.stats$orthogroup

	target.proteins = read.table('fasta/aa/GCF_003255815.1_Tgel_1.0_protein.faa.fai',sep='\t',stringsAsFactors=FALSE)[1:2]
	target.proteins = subset(target.proteins,V1 %in% sco.ortho$thegel)
	names(target.proteins) = c('tgel.protein','aa.length')

	human.proteins = read.table('fasta/aa/GCF_000001405.39_GRCh38.p13_protein.faa.fai',sep='\t',stringsAsFactors=FALSE)[1:2]
	human.proteins = subset(human.proteins,V1 %in% sco.ortho$thegel)
	names(target.proteins) = c('tgel.protein','aa.length')

	sco.ortho.tgel = subset(sco.ortho,thegel %in% target.proteins[[1]])
	rownames(target.proteins) = rownames(sco.ortho.tgel)[match(target.proteins[[1]],sco.ortho.tgel$thegel)]

	sco.stats = data.frame(sco.stats[rownames(target.proteins),],target.proteins)
	sco.stats$cds.length = with(sco.stats,aa.length * 3)
	sco.stats = sco.stats[order(rownames(sco.stats)),]

	write.table(sco.stats,file='sco_stats.txt',sep='\t',row.names=FALSE,quote=FALSE)
} else {
	sco.stats = read.delim('sco_stats.txt',stringsAsFactors=FALSE)
}

all.selection = merge(all.selection,orthogroup.dims,by='id',all.x=TRUE,all.y=FALSE)
all.selection = merge(all.selection,sco.stats,by.x='id',by.y='orthogroup',all.x=TRUE,all.y=FALSE)

all.selection$alignment.perct = with(all.selection,nchar/cds.length)




# FILTER
all.selection.backup = all.selection

i = read.table('results/paml_results_all.txt',header=FALSE,sep='\t',stringsAsFactors=FALSE)[[1]]

all.selection.strict = subset(all.selection.backup,n.duplicated < 1 & max.duplicated <= 1 & ntax >= 36 & (nchar >= 120 | alignment.perct >= 0.25))
all.selection.loose = subset(all.selection.backup,n.duplicated <= 2 & max.duplicated <= 2 & ntax >= 36 & (nchar >= 120 | alignment.perct >= 0.25))

all.selection.list = list('strict' = all.selection.strict, 'loose' = all.selection.loose)





o2go = readRDS('checkpoints/orthogroups_to_go.rds')

for (x in c('strict')) {

all.selection = all.selection.list[[x]]

all.selection.go.all = merge(all.selection,o2go,by.x='id',by.y='orthogroup',all.x=TRUE,all.y=FALSE)

# all.selection.go = subset(all.selection.go.all,namespace_1003 %in% 'biological_process')

all.selection.go = all.selection.go.all

# count.cutoff = 3
# 
# # Preserve full datasets
# go.cafe.all = go.cafe
# 
# go.cafe = subset(go.cafe,count.all >= count.cutoff)
# 
# go.cafe$overrepresentation = unlist(lapply(1:nrow(go.cafe),function(i) binom.test(go.cafe$count.expansion[i],go.cafe$count.all[i],go.cafe$prob[1],alternative='greater')$p.value))
# 
# # go.cafe$overrepresentation = rep.test(go.cafe,sum(cafe$expansion),sig.column='count.expansion')
# 
# go.cafe$overrepresentation.fdr = p.adjust(go.cafe$overrepresentation,'fdr')


# # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# #                           Overrepresentation test                             #
# # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# 
# # Binomial test. Go through each GO term and compare its count among significant genes to the expected count based on its proportional representation in the full dataset 
# rep.test = function(annotation.dataset,n.trials,alt='greater',sig.column='count.sig',prob.column='prob') {
# 	if (as.logical(n.trials)) {
# 		do.call(c,lapply(1:nrow(annotation.dataset),function(x) {
# 			binom.test(annotation.dataset[[sig.column]][x],n.trials,p=annotation.dataset[[prob.column]][x],alternative=alt)$p.value
# 		}))
# 	} else {
# 		rep(NA,nrow(annotation.dataset))
# 	}
# }

# # Enrichment test. Compare stats associated with the GO term to stats not associated with it
# enrich.test = function(annotation.dataset,gene.dataset,stat.column,alt='greater',accession.column='go.accession',annotation.column='go.annotations') {
# 	do.call(c,lapply(annotation.dataset[[accession.column]],function(x) {
# 		enrich.dataframe = data.frame(
# 			in.group = factor(gene.dataset[[annotation.column]] %in% x,levels=c('TRUE','FALSE')),
# 			stat = gene.dataset[[stat.column]]
# 		)
# 		
# 		wilcox.test(stat ~ in.group,data=enrich.dataframe,alternative=alt)$p.value
# 	}))
# }
# 
# 
# go.cafe$overrepresentation = unlist(lapply(1:nrow(go.cafe),function(i) binom.test(go.cafe$count.expansion[i],go.cafe$count.all[i],sum(cafe$expansion) / nrow(cafe),alternative="greater")$p.value))
# 
# # go.cafe$overrepresentation = rep.test(go.cafe,sum(cafe$expansion),sig.column='count.expansion')
# 
# go.cafe$overrepresentation.fdr = p.adjust(go.cafe$overrepresentation,'fdr')





library(topGO)
library(parallel)

ortho.go = all.selection.go[,c('id','go_id')]
ortho.id.to.go = mclapply(unique(ortho.go$id),function(x) sort(ortho.go$go_id[ortho.go$id == x]),mc.cores=4)
names(ortho.id.to.go) = unique(ortho.go$id)

orthogroups = all.selection$paml.p
names(orthogroups) = all.selection$id

# go.data = new('topGOdata',description='Simple session',ontology='BP',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)
go.data1 = new('topGOdata',description='Simple session',ontology='BP',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)
go.data2 = new('topGOdata',description='Simple session',ontology='MF',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)
go.data3 = new('topGOdata',description='Simple session',ontology='CC',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)

go.data = go.data1

if (!file.exists('checkpoints/orthogroup_paml_go_names.rds')) {
	ortho.topgo = data.frame(go_id = c(go.data1@graph@nodes,go.data2@graph@nodes,go.data3@graph@nodes),stringsAsFactors=FALSE)
	library(tidyverse)
	library(XML)
	
	ortho.topgo = merge(ortho.topgo,unique(subset(o2go,select=c('go_id','name_1006'))),by='go_id',all.x=TRUE)

	to.do = ortho.topgo$go_id[is.na(ortho.topgo$name_1006)]

	# Grab names and namespaces from the GO database
	library(GO.db)

	go.terms = data.frame(
		go_id = to.do,
		name_1006 = Term(to.do),
		stringsAsFactors=FALSE
	)

	ortho.topgo$name_1006[is.na(ortho.topgo$name_1006)] = go.terms[ortho.topgo$go_id[is.na(ortho.topgo$name_1006)],]$name_1006

	to.do = ortho.topgo$go_id[is.na(ortho.topgo$name_1006)]

	# Fetch missing go names from the AMIGO site
	for (i in to.do) {
		message(i)
		search.url = paste0('http://amigo.geneontology.org/amigo/term/',i)
		repeat {
			xml.common = try(search.url %>% htmlParse %>% xmlChildren %>% `[[`(3) %>% xmlChildren %>% `[[`(4) %>% xmlChildren %>% `[[`(24),silent=TRUE)
			if (!'try-error' %in% class(xml.common)) break
		}
		alternate = xml.common %>% `[[`(18) %>% `[[`(4) %>% `[[`(6) %>% `[[`(42) %>% `[[`(4) %>% xmlValue
		if (!is.na(alternate)) {
			search.url = paste0('http://amigo.geneontology.org/amigo/term/',alternate)
			xml.common = search.url %>% htmlParse %>% xmlChildren %>% `[[`(3) %>% xmlChildren %>% `[[`(4) %>% xmlChildren %>% `[[`(24)
		}
		ortho.topgo$name_1006[match(i,ortho.topgo$go_id)] = xml.common %>% `[[`(6) %>% `[[`(2) %>% `[[`('text') %>% xmlValue
	}
	saveRDS(ortho.topgo,file='checkpoints/orthogroup_paml_go_names.rds')
} else {
	ortho.topgo = readRDS('checkpoints/orthogroup_paml_go_names.rds')
}

#fetWeight01 = runTest(go.data,algorithm='weight01',statistic='fisher')
#fetParentchild = runTest(go.data,algorithm='parentchild',statistic='fisher')

# Focus on KS
ksWeight01.bp = runTest(go.data1,algorithm='weight01',statistic='ks')
#ksElim = runTest(go.data,algorithm='elim',statistic='ks')

ksWeight01.mf = runTest(go.data2,algorithm='weight01',statistic='ks')
ksWeight01.cc = runTest(go.data3,algorithm='weight01',statistic='ks')

ortho.topgo = within(ortho.topgo,{
# 	paml.fet.weight01 = fetWeight01@score[match(ortho.topgo$go_id,names(fetWeight01@score))]
# 	paml.fet.parentchild = fetParentchild@score[match(ortho.topgo$go_id,names(fetParentchild@score))]
	paml.bp.ks.weight01 = ksWeight01.bp@score[match(ortho.topgo$go_id,substr(names(ksWeight01.bp@score),1,10))]
	paml.mf.ks.weight01 = ksWeight01.mf@score[match(ortho.topgo$go_id,substr(names(ksWeight01.mf@score),1,10))]
	paml.cc.ks.weight01 = ksWeight01.cc@score[match(ortho.topgo$go_id,substr(names(ksWeight01.cc@score),1,10))]
#	paml.ks.elim = ksElim@score[match(ortho.topgo$go_id,names(ksElim@score))]
})

# ortho.topgo = within(ortho.topgo,{
# #	paml.fet.weight01.fdr = p.adjust(paml.fet.weight01,'fdr')
# #	paml.fet.parentchild.fdr = p.adjust(paml.fet.parentchild,'fdr')
# 	paml.ks.weight01.fdr = p.adjust(paml.ks.weight01,'fdr')
# 	paml.mf.ks.weight01.fdr = p.adjust(paml.mf.ks.weight01,'fdr')
# #	paml.ks.elim.fdr = p.adjust(paml.ks.elim,'fdr')
# })

orthogroups = all.selection$busted.p
names(orthogroups) = all.selection$id
orthogroups = orthogroups[!is.na(orthogroups)]

go.data1 = new('topGOdata',description='Simple session',ontology='BP',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)
go.data2 = new('topGOdata',description='Simple session',ontology='MF',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)
go.data3 = new('topGOdata',description='Simple session',ontology='CC',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)

#fetWeight01 = runTest(go.data,algorithm='weight01',statistic='fisher')
# fetParentchild = runTest(go.data,algorithm='parentchild',statistic='fisher')

# Focus on KS
ksWeight01.bp = runTest(go.data1,algorithm='weight01',statistic='ks')
# ksElim = runTest(go.data,algorithm='elim',statistic='ks')

ksWeight01.mf = runTest(go.data2,algorithm='weight01',statistic='ks')
ksWeight01.cc = runTest(go.data3,algorithm='weight01',statistic='ks')

ortho.topgo = within(ortho.topgo,{
# 	busted.fet.weight01 = fetWeight01@score[match(ortho.topgo$go_id,names(fetWeight01@score))]
# 	busted.fet.parentchild = fetParentchild@score[match(ortho.topgo$go_id,names(fetParentchild@score))]
	busted.bp.ks.weight01 = ksWeight01.bp@score[match(ortho.topgo$go_id,substr(names(ksWeight01.bp@score),1,10))]
	busted.mf.ks.weight01 = ksWeight01.mf@score[match(ortho.topgo$go_id,substr(names(ksWeight01.mf@score),1,10))]
	busted.cc.ks.weight01 = ksWeight01.cc@score[match(ortho.topgo$go_id,substr(names(ksWeight01.cc@score),1,10))]
#	busted.ks.elim = ksElim@score[match(ortho.topgo$go_id,names(ksElim@score))]
})

# ortho.topgo = within(ortho.topgo,{
# #	busted.fet.weight01.fdr = p.adjust(busted.fet.weight01,'fdr')
# #	busted.fet.parentchild.fdr = p.adjust(busted.fet.parentchild,'fdr')
# 	busted.ks.weight01.fdr = p.adjust(busted.ks.weight01,'fdr')
# 	busted.mf.ks.weight01.fdr = p.adjust(busted.mf.ks.weight01,'fdr')
# #	busted.ks.elim.fdr = p.adjust(busted.ks.elim,'fdr')
# })

orthogroups = all.selection$relax.p
names(orthogroups) = all.selection$id
orthogroups = orthogroups[!is.na(orthogroups)]

go.data1 = new('topGOdata',description='Simple session',ontology='BP',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)
go.data2 = new('topGOdata',description='Simple session',ontology='MF',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)
go.data3 = new('topGOdata',description='Simple session',ontology='CC',allGenes=orthogroups,geneSel=function(orthogroups) orthogroups<0.05,nodeSize=10,annotationFun=annFUN.gene2GO,gene2GO=ortho.id.to.go)

#fetWeight01 = runTest(go.data,algorithm='weight01',statistic='fisher')
#fetParentchild = runTest(go.data,algorithm='parentchild',statistic='fisher')

# Focus on KS
ksWeight01.bp = runTest(go.data1,algorithm='weight01',statistic='ks')
ksWeight01.mf = runTest(go.data2,algorithm='weight01',statistic='ks')
ksWeight01.cc = runTest(go.data3,algorithm='weight01',statistic='ks')
#ksElim = runTest(go.data,algorithm='elim',statistic='ks')

ortho.topgo = within(ortho.topgo,{
 #	relax.fet.weight01 = fetWeight01@score[match(ortho.topgo$go_id,names(fetWeight01@score))]
 #	relax.fet.parentchild = fetParentchild@score[match(ortho.topgo$go_id,names(fetParentchild@score))]
	relax.bp.ks.weight01 = ksWeight01.bp@score[match(ortho.topgo$go_id,substr(names(ksWeight01.bp@score),1,10))]
	relax.mf.ks.weight01 = ksWeight01.mf@score[match(ortho.topgo$go_id,substr(names(ksWeight01.mf@score),1,10))]
	relax.cc.ks.weight01 = ksWeight01.cc@score[match(ortho.topgo$go_id,substr(names(ksWeight01.cc@score),1,10))]
#	relax.ks.elim = ksElim@score[match(ortho.topgo$go_id,names(ksElim@score))]
})

# ortho.topgo = within(ortho.topgo,{
# #	relax.fet.weight01.fdr = p.adjust(relax.fet.weight01,'fdr')
# #	relax.fet.parentchild.fdr = p.adjust(relax.fet.parentchild,'fdr')
# 	relax.ks.weight01.fdr = p.adjust(relax.ks.weight01,'fdr')
# #	relax.ks.elim.fdr = p.adjust(relax.ks.elim,'fdr')
# })

ortho.topgo$namespace_1003 = Ontology(ortho.topgo$go_id)

ortho.topgo = within(ortho.topgo,{
	paml.bp.ks.weight01.fdr    = p.adjust(paml.bp.ks.weight01  ,'fdr')
	paml.mf.ks.weight01.fdr    = p.adjust(paml.mf.ks.weight01  ,'fdr')
	paml.cc.ks.weight01.fdr    = p.adjust(paml.cc.ks.weight01  ,'fdr')
	busted.bp.ks.weight01.fdr  = p.adjust(busted.bp.ks.weight01,'fdr')
	busted.mf.ks.weight01.fdr  = p.adjust(busted.mf.ks.weight01,'fdr')
	busted.cc.ks.weight01.fdr  = p.adjust(busted.cc.ks.weight01,'fdr')
	relax.bp.ks.weight01.fdr   = p.adjust(relax.bp.ks.weight01 ,'fdr')
	relax.mf.ks.weight01.fdr   = p.adjust(relax.mf.ks.weight01 ,'fdr')
	relax.cc.ks.weight01.fdr   = p.adjust(relax.cc.ks.weight01 ,'fdr')
})

ortho.topgo = ortho.topgo[order(ortho.topgo$namespace_1003,ortho.topgo$paml.bp.ks.weight01,ortho.topgo$paml.mf.ks.weight01,ortho.topgo$paml.cc.ks.weight01),]
rownames(ortho.topgo) = NULL
saveRDS(ortho.topgo,file=paste0('checkpoints/orthogroup_togo_results_',x,'.rds'))


o2g = readRDS('checkpoints/orthogroups_to_genes.rds')
o2g.sco = do.call(cbind,mclapply(o2g,function(x) {
	unlist(lapply(x[all.selection$id],function(y) {
		if (length(y) > 1) NA else y		
	}))
},mc.cores=4))

library(biomaRt)

o2g.sco.target = o2g.sco[,'thegel']
tgel = useEnsembl(biomart='ENSEMBL_MART_ENSEMBL',dataset='tgelada_gene_ensembl',version=102)
gene.info = getBM(attributes=c('ensembl_gene_id','external_gene_name','entrezgene_id'),filters='entrezgene_id',values=o2g.sco.target,mart=tgel)

gene.info$gene = gene.info$external_gene_name
gene.info$gene[is.na(gene.info$gene) | !nchar(gene.info$gene)] = gene.info$ensembl_gene_id[is.na(gene.info$gene) | !nchar(gene.info$gene)]
gene.info$id = names(o2g.sco.target[match(gene.info$entrezgene_id,o2g.sco.target)])

gene.name = subset(gene.info,select=c('id','gene'))

o2g.sco.target = o2g.sco[,'homsap']
hsap = useEnsembl(biomart='ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl',version=102)
gene.info = getBM(attributes=c('ensembl_gene_id','external_gene_name','entrezgene_id'),filters='entrezgene_id',values=o2g.sco.target,mart=hsap)

gene.info$gene = gene.info$external_gene_name
gene.info$gene[is.na(gene.info$gene) | !nchar(gene.info$gene)] = gene.info$ensembl_gene_id[is.na(gene.info$gene) | !nchar(gene.info$gene)]
gene.info$id = names(o2g.sco.target[match(gene.info$entrezgene_id,o2g.sco.target)])

gene.name = rbind(gene.name,subset(gene.info,!id %in% gene.name$id,select=c('id','gene')))

gene.name = do.call(rbind,lapply(split(gene.name,gene.name$id),function(x) {
	if (any(!grepl('^ENS',x$gene) & nchar(x$gene) != 18)) x = subset(x,!grepl('^ENS',gene) & nchar(gene) != 18)
	gene = paste(unique(x$gene),collapse=',')
	x$gene = gene
	unique(x)
}))

all.selection = merge(all.selection,gene.name,by='id',all.x=TRUE)

all.selection$paml.q = p.adjust(all.selection$paml.p,'fdr')
all.selection$busted.q = p.adjust(all.selection$busted.p,'fdr')

all.selection = all.selection[order(all.selection$id),]
rownames(all.selection) = NULL

# subset(ortho.topgo,paml.ks.weight01.fdr < 0.05 & busted.ks.weight01.fdr < 0.05)
# subset(all.selection,id %in% subset(subset(all.selection.go,go_id %in% 'GO:0001666'),paml.p < 0.05)$id)# response to hypoxia
# subset(all.selection,id %in% subset(subset(all.selection.go,go_id %in% 'GO:0071456'),paml.p < 0.05)$id)# cellular response to hypoxia
# subset(all.selection,id %in% subset(subset(all.selection.go,go_id %in% 'GO:0001525'),paml.p < 0.05)$id)# angiogenesis
# subset(all.selection,id %in% subset(subset(all.selection.go,go_id %in% 'GO:0001701'),paml.p < 0.05)$id)# in utero embryonic development
# subset(all.selection,id %in% subset(subset(all.selection.go,go_id %in% 'GO:0006979'),paml.p < 0.05)$id)# response to oxidative stress
# subset(all.selection,id %in% subset(subset(all.selection.go,go_id %in% 'GO:0007626'),paml.p < 0.05)$id)# locomotory behavior


saveRDS(all.selection,file=paste0('checkpoints/selection_sco_results_',x,'.rds'))

}

# all.selection.strict.human = subset(all.selection.strict,!is.na(hsap.protein))
library(biomaRt)

hsap = useEnsembl(biomart='ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl',version=102)

entrez.human = unlist(o2g$homsap)
all.selection.strict$entrezgene_id = entrez.human[all.selection.strict$id]

subset(all.selection.strict,p.adjust(paml.p,'fdr') < 0.2)

gene.info = getBM(attributes=c('ensembl_gene_id','entrezgene_id','hgnc_id'),filters='entrezgene_id',values=all.selection.strict$entrezgene_id,mart=hsap)

all.selection.strict = merge(all.selection.strict,gene.info,by='entrezgene_id',all.x=TRUE)

subset(all.selection.strict,!is.na(ensembl_gene_id) & p.adjust(paml.p,'fdr') < 0.2)$hgnc_id