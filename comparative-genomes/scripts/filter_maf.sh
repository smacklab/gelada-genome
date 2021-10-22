#!/bin/bash

module load maffilter/1.3.1

 maffilter input.file=test.maf input.file.compression=none 'maf.filter=WindowSplit(preferred_size=100,align=adjust,keep_small_blocks=yes),SequenceStatistics(statistics=(BlockLength,AlnScore,BlockCounts,SiteStatistics(species=(homo_sapiens,canis_familiaris,theropithecus_gelada))),file=test.csv,ref_species=homo_sapiens)'

maffilter \
	input.file=test.maf \
	input.file.compression=gzip \
	
maf.filter=                                 \
    AlnFilter2(                             \
        species=(species1,species2,etc),    \
        window.size=10,                     \
        window.step=1,                      \
        max.gap=1,                          \
        max.pos=1,                          \
        relative=no,                        \
        missing_as_gap=yes,                 \
        file=data.trash_aln.maf.gz,         \
        compression=gzip)                  


maffilter input.file=46_mammals.epo.2_24.maf.gz input.file.compression=gzip maf.filter='AlnFilter2(species=(macaca_mulatta,homo_sapiens),window.size=1,window.step=1,max.gap=0,file=trash.maf.gz,compression=gzip),Output(file=test.filtered.maf.gz,compression=gzip,mask=yes)'