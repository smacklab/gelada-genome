#!/usr/bin/env Rscript

arguments = commandArgs(trailing=TRUE)

window.size = as.integer(arguments[1])
percent.intact = as.integer(arguments[2])
percent.identity = as.integer(arguments[3])

if (length(arguments) > 3) {
	n.cores = as.integer(arguments[4])
} else {
	n.cores = 24
}

source('config/settings.cfg')

cons = readRDS(file.path('checkpoints',paste0('accelerated_regions_cons_',maf_version,'.window',window.size,'.intact',percent.intact,'.identity',percent.identity,'.rds')))

stats = readRDS(file.path('checkpoints',paste0('accelerated_regions_stats_',maf_version,'.window',window.size,'.intact',percent.intact,'.identity',percent.identity,'.rds')))

cons.all = merge(cons,stats,by='maf')

cons.pass = subset(cons.all,SequenceLengthForhomo_sapiens >= 45 & SequenceLengthForrattus_norvegicus >= 45 & SequenceLengthFormus_musculus >= 45)

cons.sig = subset(cons.pass,qval < 0.2 & SequenceLengthFortheropithecus_gelada >= 1)

gff = read.delim(gzfile('data/Theropithecus_gelada.Tgel_1.0.104.gff3.gz'),comment.char='#',header=FALSE)
gff = subset(gff,V1 %in% cons.all$chr)
library(GenomicRanges)

i.cons = with(cons.sig,GRanges(chr,IRanges(Start,Stop,names=maf),'*'))
i.gff = with(gff,GRanges(V1,IRanges(V4,V5),'*'))

gff.olaps = findOverlaps(i.cons,i.gff)

gff.hits = data.frame(cons.sig[queryHits(gff.olaps),],
	gff[subjectHits(gff.olaps),c('V1','V3','V4','V5')])

# Fraction in protein-coding gene vs. no
table(unlist(lapply(split(gff.hits,gff.hits$maf),function(x) if ('CDS' %in% x$V3) 'exonic: CDS' else if ('five_prime_UTR' %in% x$V3) 'exonic: five_prime_UTR' else if ('three_primer_UTR' %in% x$V3) 'exonic: three_primer_UTR' else if ('gene' %in% x$V3) 'intronic' else 'intergenic')))

type.hits = reshape2::melt(unlist(lapply(split(gff.hits,gff.hits$maf),function(x) if ('CDS' %in% x$V3) 'exonic: CDS' else if ('five_prime_UTR' %in% x$V3) 'exonic: five_prime_UTR' else if ('three_primer_UTR' %in% x$V3) 'exonic: three_primer_UTR' else if ('gene' %in% x$V3) 'intronic' else 'intergenic')))




# Find nearest gene
gff.gene = subset(gff,V3 == 'gene')

i.gene = with(gff.gene,GRanges(V1,IRanges(V4,V5),'*'))

gene.olaps = nearest(i.cons,i.gene)

gene.hits = data.frame(cons.sig,
	gff.gene[gene.olaps,c('V1','V3','V4','V5','V9')])

gene.hits$ensembl_gene_id = gsub('.+?gene:(ENSTGEG[0-9]{11}).*','\\1',gene.hits$V9)



# Now match all cons.pass to nearest genes
i.pass = with(cons.pass,GRanges(chr,IRanges(Start,Stop,names=maf),'*'))
pass.olaps = nearest(i.pass,i.gene)
pass.hits = data.frame(cons.pass,
	gff.gene[pass.olaps,c('V1','V3','V4','V5','V9')])
pass.hits$ensembl_gene_id = gsub('.+?gene:(ENSTGEG[0-9]{11}).*','\\1',pass.hits$V9)



library(biomaRt)
tgel = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',
	host='archive.ensembl.org',
	version=101,
	dataset='tgelada_gene_ensembl')

tgel.go = getBM(
	attributes=c('ensembl_gene_id','go_id'),
	mart = tgel)

tgel.go = subset(tgel.go,!is.na(go_id) & nchar(go_id))

gene2go = split(tgel.go$go_id,tgel.go$ensembl_gene_id)

library(topGO)

all.genes = integer(length(unique(pass.hits$ensembl_gene_id)))
names(all.genes) = unique(pass.hits$ensembl_gene_id)

all.genes[names(all.genes) %in% gene.hits$ensembl_gene_id] = 1

all.fet.topgo = new('topGOdata',
	description='Simple session',
	ontology='BP',
	allGenes=all.genes,
	geneSelectionFun=function(x) x > 0,
	nodeSize = 10,
	annotationFun = annFUN.gene2GO,
	gene2GO = gene2go
)

all.fet.test = runTest(all.fet.topgo,algorithm='parentchild',statistic='fisher')

# Grab names and namespaces from the GO database
library(GO.db)

go.terms = data.frame(
	go_id = names(all.fet.test@score),
	go_namespace = Ontology(names(all.fet.test@score)),
	go_name = Term(names(all.fet.test@score)),
	stringsAsFactors=FALSE
)
rownames(go.terms) = go.terms$go_id

all.fet.results = reshape2::melt(all.fet.test@score)
all.fet.results$go_id = rownames(all.fet.results)

all.fet.results = merge(all.fet.results,go.terms,by='go_id',all.x=TRUE)
all.fet.results = all.fet.results[order(all.fet.results$value),]
names(all.fet.results)[names(all.fet.results) == 'value'] = 'pval'
all.fet.results$qval = p.adjust(all.fet.results$pval,'fdr')




# Now greatly relax the criteria for GARs and try again
all.genes = integer(length(unique(pass.hits$ensembl_gene_id)))
names(all.genes) = unique(pass.hits$ensembl_gene_id)

all.genes[names(all.genes) %in% subset(pass.hits,pval < 0.01)$ensembl_gene_id] = 1

loose.fet.topgo = new('topGOdata',
	description='Simple session',
	ontology='BP',
	allGenes=all.genes,
	geneSelectionFun=function(x) x > 0,
	nodeSize = 10,
	annotationFun = annFUN.gene2GO,
	gene2GO = gene2go
)

loose.fet.test = runTest(loose.fet.topgo,algorithm='parentchild',statistic='fisher')

loose.fet.results = reshape2::melt(loose.fet.test@score)
loose.fet.results$go_id = rownames(loose.fet.results)

loose.fet.results = merge(loose.fet.results,go.terms,by='go_id',all.x=TRUE)
loose.fet.results = loose.fet.results[order(loose.fet.results$value),]
names(loose.fet.results)[names(loose.fet.results) == 'value'] = 'pval'
loose.fet.results$qval = p.adjust(loose.fet.results$pval,'fdr')



# Simply match GO BPs to significant GARs
sig.hits = subset(pass.hits,maf %in% cons.sig$maf)

go.sig = data.frame(
	go_id = unique(na.omit(unlist(gene2go[unique(sig.hits$ensembl_gene_id)]))),
	go_namespace = Ontology(unique(na.omit(unlist(gene2go[unique(sig.hits$ensembl_gene_id)])))),
	go_name = Term(unique(na.omit(unlist(gene2go[unique(sig.hits$ensembl_gene_id)])))),
	stringsAsFactors=FALSE
)

sig.go = do.call(rbind,lapply(split(sig.hits,1:nrow(sig.hits)),function(x) {
	go.ids = gene2go[[x$ensembl_gene_id]]
	if (is.null(go.ids)) {
		go.string = ''
	} else {
		go.string = paste(subset(go.sig,go_id %in% go.ids & go_namespace == 'BP')$go_name,collapse='|')
	}
	data.frame(x,go_matches=go.string)
}))





cons.sig = cons.sig[order(as.integer(gsub('[a-z]$','',cons.sig$chr)),gsub('[0-9]','',cons.sig$chr),cons.sig$pos),]
cons.sig$gar_id = paste0('GAR',1:nrow(cons.sig))


cons.merge = subset(cons.sig,select=c('gar_id','block.x'))
names(cons.merge)[2] = 'block'

phylop = readRDS(file.path('checkpoints',paste0('accelerated_regions_phylop_',maf_version,'.window',window.size,'.intact',percent.intact,'.identity',percent.identity,'.rds')))

cons.sig.phylop = subset(phylop,block %in% cons.sig$block.x)

cons.scores = cons.sig$pval
names(cons.scores) = cons.sig$block.x
cons.sig.phylop$cons.pval = cons.scores[cons.sig.phylop$block]

cons.sig.phylop = within(cons.sig.phylop,{
	acc.pval = pval
	acc.pval[pval>0]=-1;
	acc.pval=-acc.pval
	block=factor(block,levels=cons.sig[order(cons.sig$pval),]$block.x)
})

cons.sig.phylop = merge(cons.sig.phylop,cons.merge,by='block',all.x=TRUE)

cons.sig.phylop$gar_id = factor(cons.sig.phylop$gar_id,levels=paste0('GAR',1:nrow(cons.sig)))

library(ggplot2)
p = ggplot(cons.sig.phylop,aes(pos,-log10(acc.pval))) +
	geom_point(size=0.01,shape=19) +
	facet_wrap(~gar_id,nrow=3,scales='free_x') +
	theme_classic(base_size=16) +
	theme(
		strip.background=element_blank(),
#		strip.text=element_blank(),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.line.x=element_blank()
	) +
	xlab('Position') +
	ylab(expression(-log[10]~italic(p)))
ggsave(p,file='gar_phylop.pdf',useDingbats=FALSE,height=5,width=9)

# Bring in ChromHMM states
chrom.hmm = readRDS('checkpoints/GAR_chromHMMstates.rds')
# Br

chrom.hmm.set = subset(chrom.hmm,
	tissuetype %in% c('Lung','Spleen','Liver','Primary hematopoietic stem cells','Osteoblast Primary Cells','Ovary','Placenta','Placenta Amnion','Esophagus','Aorta','Gastric','Thymus') |
	grepl('^Brain',tissuetype) |
	grepl('Atrium$',tissuetype) |
	grepl('Ventricle$',tissuetype) |
	grepl('^Fetal',tissuetype) |
	grepl('^Stomach',tissuetype) |
	grepl('Intestine',tissuetype) |
	grepl('Duodenum',tissuetype) |
	grepl('Colon',tissuetype) |
	grepl('^Pancrea',tissuetype) |
	grepl('^Skeletal Muscle',tissuetype) |
	grepl('Muscle$',tissuetype) |
	grepl('Primary Cells$',tissuetype) |
	grepl('^Primary.+?blood.*',tissuetype)
)
# Check types that were not picked
# setdiff(unique(chrom.hmm$tissuetype),chrom.hmm.set$tissuetype)


# chrom.hmm.set = subset(chrom.hmm,TRUE)

# chrom.hmm.marks = as.data.frame(subset(chrom.hmm.set,state %in% c('Enhancers','Strong transcription','Genic enhancers','Active TSS','Flanking Active TSS','Bivalent/Poised TSS','Bivalent Enhancer')))

regulatory.states = c(
	'Active TSS',
	'Flanking Active TSS',
	'Transcr. at gene 5\' and 3\'',
#	'Strong transcription',
#	'Weak transcription',
	'Genic enhancers',
	'Enhancers',
#	'ZNF genes & repeats',
#	'Heterochromatin',
	'Bivalent/Poised TSS',
	'Flanking Bivalent TSS/Enh',
	'Bivalent Enhancer'#,
#	'Repressed PolyComb',
#	'Weak Repressed PolyComb',
#	'Quiescent/Low'
)

chrom.hmm.marks = as.data.frame(subset(chrom.hmm.set,state %in% regulatory.states))

# Summarize by GAR
chrom.hmm.summary = do.call(rbind,lapply(split(chrom.hmm.marks,chrom.hmm.marks$GAR),function(x) {
	chrom.hmm.string = paste(unlist(lapply(split(x,x$state),function(y) {
		paste0(unique(y$state),' (',paste(y$tissuetype,collapse='|'),')')
	})),collapse='; ')
	data.frame(maf=unique(x$GAR),chrom_hmm=chrom.hmm.string)
}))

sig.summary = merge(sig.go,chrom.hmm.summary,by='maf',all.x=TRUE)

sig.summary = sig.summary[with(sig.summary,order(as.integer(gsub('[A-z]*','',chr)),Start)),]
sig.summary = sig.summary[with(sig.summary,order(pval)),]

sig.summary$element_type = type.hits[sig.summary$maf,'value']

# Add gene names

tgel.names = getBM(
	attributes=c('ensembl_gene_id','external_gene_name'),
	mart = tgel)

sig.name = subset(tgel.names,ensembl_gene_id %in% sig.summary$ensembl_gene_id)
name.key = sig.name$external_gene_name
names(name.key) = sig.name$ensembl_gene_id

sig.summary$external_gene_name = name.key[sig.summary$ensembl_gene_id]



sig.summary = subset(sig.summary,select=c(
	'maf','str','lnl_null','lnl_alt','pval','qval','Chr','Start','Stop','BlockSize','BlockLength','SequenceLengthFortheropithecus_gelada','SequenceLengthForhomo_sapiens','homo_sapiens.Chr','homo_sapiens.Start','homo_sapiens.Stop','ensembl_gene_id','go_matches','chrom_hmm','element_type','external_gene_name'
))

dir.create('tables',showWarnings=FALSE)

write.table(sig.summary,file='tables/gars_acceleration.txt',sep='\t',row.names=FALSE,quote=FALSE)

sig.summary = sig.summary[with(sig.summary,order(as.integer(gsub('[A-z]*','',Chr)),gsub('[0-9]*','',Chr),Start)),]
write.table(sig.summary,file='tables/gars_position.txt',sep='\t',row.names=FALSE,quote=FALSE)
