#!/usr/bin/env Rscript

# Make protein to taxon join list, but this isn't actually needed
all.proteins = read.table('id_lists/taxon_joins/prot2taxon.txt',sep='\t',header=FALSE,stringsAsFactors=FALSE)
saveRDS(all.proteins,file='id_lists/taxon_joins/prot2taxon.rds')

all.assemblies = read.table('assembly_list.txt',sep='\t',stringsAsFactors=FALSE)

assemblies = all.assemblies$V4
names(assemblies) = gsub('-','.',gsub(' ','_',with(all.assemblies,paste0(V1,'_',V2,'_protein'))))

reference.assembly = with(all.assemblies[1,],paste0(V1,'_',V2))
reference.assembly.short.name = all.assemblies$V4[1]
p2g = read.table(paste0('id_lists/feature_joins/',reference.assembly,'_prot2gene.txt'),sep='\t')

orthogroups = read.delim('id_lists/orthogroups/orthogroups.tsv',stringsAsFactors=FALSE)
rownames(orthogroups) = orthogroups$Orthogroup
orthogroups$Orthogroup = NULL

names(orthogroups) = assemblies[names(orthogroups)]

library(stringr)

orthogroup.sizes = do.call(data.frame,lapply(orthogroups,function(x) ifelse(nchar(x),str_count(x,', ')+1,0)))
rownames(orthogroup.sizes) = NULL

orthogroup.sizes = data.frame(
	Description = rownames(orthogroups),
	ID = rownames(orthogroups),
	orthogroup.sizes
)

write.table(
	orthogroup.sizes,
	file='id_lists/orthogroups/protein_counts.txt',
	sep='\t',
	row.names=FALSE,
	quote=FALSE)
