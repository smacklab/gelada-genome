#!/usr/bin/env Rscript

arguments = commandArgs(trailingOnly = TRUE)

accession = arguments[1]
id = arguments[2]

if (length(arguments) > 2 ) db = arguments[3] else db = 'unknown'

assembly = paste(accession,id,sep='_')

fai = read.table(
		paste0('fasta/aa/',assembly,'_protein.faa.fai'),
		sep='\t',
		header=FALSE,
		stringsAsFactors=FALSE)[1:2]

names(fai) = c('protein.accession','protein.length')

if (db == 'ensembl') {

	latin.name = unlist(strsplit(accession,'_'))
	ensembl.code = paste0(substr(tolower(latin.name[1]),1,1),latin.name[length(latin.name)])

	library(biomaRt)
	mart = useMart(
			biomart='ENSEMBL_MART_ENSEMBL',
			dataset=paste0(ensembl.code,'_gene_ensembl'))
	p2all = getBM(
			attributes=c('ensembl_gene_id','ensembl_peptide_id_version','ensembl_transcript_id_version'),
			filters='ensembl_peptide_id_version',
			values=fai$protein.accession,mart=mart)

	p2g = unique(p2all[c('ensembl_gene_id','ensembl_peptide_id_version')])
	names(p2g) = c('gene.id','protein.accession')
	write.table(
			p2g,
			file=paste0('id_lists/feature_joins/',assembly,'_prot2gene.txt'),
			sep='\t',
			row.names=FALSE,
			quote=FALSE,
			col.names=FALSE)

	p2c = unique(p2all[c('ensembl_peptide_id_version','ensembl_transcript_id_version')])
	names(p2c) = c('protein.accession','cds.id')
	write.table(
			p2c,
			file=paste0('id_lists/feature_joins/',assembly,'_prot2cds.txt'),
			sep='\t',
			row.names=FALSE,
			quote=FALSE,
			col.names=FALSE)

} else if (db == 'refseq') {

	p2g = read.table(
			paste0('id_lists/feature_joins/',assembly,'_prot2gene.txt'),
			sep='\t',
			header=FALSE,
			stringsAsFactors=FALSE,
			col.names=c('gene.id','protein.accession'))

}

p2g = merge(p2g,fai,by='protein.accession',all.x=TRUE)

p2g.split = split(p2g,p2g$gene.id)

longest.isoform = do.call(rbind,lapply(p2g.split,function(x) {
	x[which.max(x$protein.length),]
}))

write(
		longest.isoform$protein.accession,
		sep='\n',
		file=paste0('id_lists/aa_longest_isoform/',assembly,'_protein.faa.txt'))
