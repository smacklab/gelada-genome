#!/usr/bin/env Rscript

library(parallel)

# Read in assembly metadata
assembly.list = read.table('assembly_list.txt',sep='\t',stringsAsFactors=FALSE)
assembly.list$assembly = gsub(' ','_',with(assembly.list,paste(V1,V2,sep='_')))

assembly.codes = assembly.list$V4
names(assembly.codes) = gsub(' ','_',assembly.list$V2)

# List files with best HMM matches per protein per assembly
best.matches = list.files(path='treefam_tools/treefam_scan/results',pattern='*_bestfams.txt')

# Read each table and add the assembly name as a new column
best.fams = do.call(rbind,mclapply(best.matches,function(x) {
	out = read.table(file.path('treefam_tools/treefam_scan/results',x),header=TRUE,sep=' ',stringsAsFactors=FALSE)
	out$assembly = gsub('treefam_hmm_(.+?)_bestfams.txt','\\1',x)
	out
},mc.cores=16))

# Use short codes for assembly columns
best.fams$assembly = as.character(assembly.codes[best.fams$assembly])

# Factorize the assembly column (so that tabular counts will report all levels from data subsets)
best.fams$assembly = factor(best.fams$assembly)

# Rename column headings
names(best.fams) = c('refseq_protein','treefam','bit_score','assembly')

# Order by assembly and protein ID
best.fams = best.fams[order(best.fams$assembly,best.fams$refseq_protein),]

# Split dataset by treefam
best.treefam = split(best.fams,best.fams$treefam)

# For each treefam, calculate and return gene/protein counts per assembly
best.treefam = do.call(rbind,mclapply(best.treefam,function(x) {
	as.data.frame.matrix(t(table(x$assembly)))
},mc.cores=16))

# Save treefam IDs
best.treefam$treefam = rownames(best.treefam)
best.treefam$notes = rownames(best.treefam)

# Reorder columns
best.treefam = best.treefam[,c('treefam','notes',names(best.treefam)[1:(ncol(best.treefam)-2)])]
names(best.treefam)[1:2] = c('Description','ID')
names(best.treefam)[3:ncol(best.treefam)] = gsub('_','',tolower(names(best.treefam)[3:ncol(best.treefam)]))
rownames(best.treefam) = NULL

write.table(best.treefam,file='id_lists/orthogroups/treefam_counts.txt',row.names=FALSE,quote=FALSE,sep='\t')
