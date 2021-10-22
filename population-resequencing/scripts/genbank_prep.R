#!/usr/bin/env Rscript

mt.genomes = scan('mt_genomes.fa',what='',sep='\n')

mt.genomes = data.frame(
	id=gsub('^>','',mt.genomes[seq(1,length(mt.genomes),2)]),
	fasta=mt.genomes[seq(2,length(mt.genomes),2)]
)

mt.genomes.split = split(mt.genomes,mt.genomes$fasta)
names(mt.genomes.split) = NULL

mt.ids = do.call(rbind,lapply(mt.genomes.split,function(x) {
	out = data.frame(id=paste(x$id,collapse=','),fasta=unique(x$fasta))
}))

rownames(mt.ids) = NULL

rownames(mt.ids) = mt.ids$id

mt.ids = mt.ids[with(mt.ids,order(-grepl('FJ785426.1',id),-grepl('^FIL',id),-grepl('^GUA',id),-grepl('^FRZ',id),id)),]

mt.ids$new.id = c('FJ785426.1',paste0(rep('hamadryas',sum(grepl('^FIL',mt.ids$id))),formatC(1:sum(grepl('^FIL',mt.ids$id)),width=2,flag=0)),paste0('gelada',formatC(1:(nrow(mt.ids)-sum(grepl('^FIL',mt.ids$id))-1),width=2,flag=0)))

translation.key = mt.ids[c('new.id','id')]

write.table(translation.key,col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t',file='genbank_seq_membership.txt')

# Write fasta

fasta.out = as.matrix(mt.ids[c('new.id','fasta')])
dimnames(fasta.out) = NULL
fasta.out[,1] = paste0('>',fasta.out[,1])

fasta.out = as.character(t(fasta.out))

write(fasta.out,file='mt_genbank_unaligned.fa',sep='\n')

# Now write one sequence per file

dir.create('mt_genbank',showWarnings=FALSE)

for (i in 1:(length(fasta.out)/2)) {
	j = 2 * i - 1
	write(fasta.out[j:(j+1)],file=file.path('mt_genbank',paste0(gsub('>','',fasta.out[j]),'.fa')),sep='\n')
}