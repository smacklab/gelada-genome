#!/usr/bin/env Rscript

ignore.checkpoints = FALSE

if (!file.exists('data/PTHR15.0_human')) {
	download.file('ftp://ftp.pantherdb.org/sequence_classifications/15.0/PANTHER_Sequence_Classification_files/PTHR15.0_human',destfile='data/PTHR15.0_human')
}

if (!file.exists('data/panther_15.0_human_hif.txt')) {
	system('cut -f 1,10 data/PTHR15.0_human | grep P00030 > data/panther_15.0_human_hif.txt')
	#system('cut -f 1,10 data/PTHR15.0_macacque > data/panther_15.0_macaque_all.txt')
}

hif = read.table('data/panther_15.0_human_hif.txt',sep='\t',comment.char='',stringsAsFactors=FALSE)

hif$hgnc = NA
hif$hgnc[grep('HGNC=',hif$V1)] = gsub('^.+?HGNC=([0-9]+).*','\\1',hif$V1)[grep('HGNC=',hif$V1)]
hif$ensembl = NA
hif$ensembl[grep('Ensembl=',hif$V1)] = gsub('^.+?Ensembl=(ENSG[0-9]+).*','\\1',hif$V1)[grep('Ensembl=',hif$V1)]

# If any ensembl IDs are missing, find them via a link to the gene name
if (ignore.checkpoints || !file.exists('checkpoints/panther_hsap_biomart.rds')) {
	library(biomaRt)
	hsap = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='hsapiens_gene_ensembl')
	hsap.genes = getBM(
		attributes=c('ensembl_gene_id','external_gene_name','hgnc_id','hgnc_symbol','chromosome_name'),
		mart = hsap)
	saveRDS(hsap.genes,file='checkpoints/panther_hsap_biomart.rds')
} else {
	message('Checkpoint found!\nLoading gene annotations from file.')
	hsap.genes = readRDS('checkpoints/panther_hsap_biomart.rds')
}

# Subset to main chromosomes to avoid haplotyes
hsap.hif = subset(hsap.genes,chromosome_name %in% 1:22 & (hgnc_id %in% with(hif,paste0('HGNC:',hgnc[!is.na(hgnc)])) | ensembl_gene_id %in% with(hif,ensembl[!is.na(ensembl)])))

digestive.enyzmes = read.table('data/Homo_protein_ids_digestive_enzymes.txt',sep='\t',stringsAsFactors=FALSE)

if (ignore.checkpoints || !file.exists('checkpoints/digestive_hsap_biomart.rds')) {
	library(biomaRt)
	hsap = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='hsapiens_gene_ensembl')
	hsap.genes = getBM(
		attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'),
		mart = hsap)
	saveRDS(hsap.genes,file='checkpoints/digestive_hsap_biomart.rds')
} else {
	message('Checkpoint found!\nLoading gene annotations from file.')
	hsap.genes = readRDS('checkpoints/digestive_hsap_biomart.rds')
}

hsap.digestive = subset(hsap.genes,external_gene_name %in% digestive.enyzmes$V2 & chromosome_name %in% 1:22)

if (!all(digestive.enyzmes$V2 %in% hsap.digestive$external_gene_name)) warning('Not all digestive genes found')

# Theropithecus (Tgel_1.0)
# Papio (Panubis1.0)
# Rhesus (Mmul_10)
# Pan (Clint_PTRv2)
# Homo (GRCh38.p13)
# Mus (GRCm38.p6)
# Goat (ARS1)
# Cow (ARS-UCD1.2)
# Cat (Felis_catus_9.0)
# Bat (mRhiFer1_v1.p)

ortholog.species = c('tgelada','panubis','mmulatta','ptroglodytes','mmusculus','chircus','btaurus','fcatus','rferrumequinum')

hsap.allgenes = c(hsap.digestive$ensembl_gene_id,hsap.hif$ensembl_gene_id)

if (ignore.checkpoints || !file.exists('checkpoints/hsap_orthologs_biomart.rds')) {
	library(biomaRt)
	hsap = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',
		dataset='hsapiens_gene_ensembl')
	hsap.orthologs = getBM(
		attributes=c(
			'ensembl_gene_id',
			'external_gene_name',
			paste0(ortholog.species,'_homolog_ensembl_gene'),
			paste0(ortholog.species,'_homolog_orthology_type')
		),
		filters = 'ensembl_gene_id',
		values = hsap.allgenes,
		mart = hsap)
	saveRDS(hsap.orthologs,file='checkpoints/hsap_orthologs_biomart.rds')
} else {
	message('Checkpoint found!\nLoading gene orthologs from file.')
	hsap.orthologs = readRDS('checkpoints/hsap_orthologs_biomart.rds')
}
