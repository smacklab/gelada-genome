#!/usr/bin/env Rscript

i = commandArgs(trailingOnly=TRUE)[1]

source('scripts/custom_functions.R')

assemblies = system('cut -f 4 assembly_list.txt',intern=TRUE)

tree = scan('data/tree.nwk',what='',quiet=TRUE)

include = gsub('>','',system(paste0('grep ">" fasta/sco/cds_corrected/',i,'.fna'),intern=TRUE))

exclude = setdiff(assemblies,include)

tree.out = tree

if (length(exclude)) {
	for (j in exclude) {
		tree.out = prune.branch(j,tree.out)
	}
}

tree.out = gsub('(thegel:[0-9.]+)','\\1{Foreground}',tree.out)

system('mkdir -p hyphy/trees')

write(tree.out,file=paste0('hyphy/trees/',i,'.nwk'))
