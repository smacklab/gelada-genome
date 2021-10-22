#!/usr/bin/env Rscript

fai = read.table('fasta/all_species/cds_corrected.fna.fai',sep='\t')

# The sixth column is the position of the fasta header
fai$V6 = as.integer(system('grep -n ">" fasta/all_species/cds_corrected.fna | cut -d ":" -f 1',intern=TRUE))

# The seventh column is the line number of the final line
fai$V7 = with(fai,V6 + ceiling(V2/V4))

# Identify the rows that are not multiples of 3
fai$V8 = fai$V2 %% 3

# Last codon has one nucleotide
fai.1 = subset(fai,V8 == 1)

# Last codon has two nucleotides
fai.2 = subset(fai,V8 == 2)

# If the final codon has two nucleotides, tack on one ambiguous base (N) to the end
if (nrow(fai.2)) {
	for (i in 1:nrow(fai.2)) {
		message('Fixing CDS sequence ',fai.2$V1[i],' by tacking on one ambiguous base.')
		# Take into account the case in which adding a base exceeds the base number per line
# 		if (nchar(system(paste0('sed -n ',fai.2$V7[i],'p fasta/all_species/cds_corrected_all_codons.fna'),intern=TRUE)) >= fai.2$V4[i]) {
# 			system(paste0('sed -i "',fai.2$V7[i]+1,' i N" fasta/all_species/cds_corrected_all_codons.fna'))
# 			fai.2$V7[fai.2$V7 > fai.2$V7[i]] = fai.2$V7[fai.2$V7 > fai.2$V7[i]] + 1
# 			fai.1$V7[fai.1$V7 > fai.2$V7[i]] = fai.1$V7[fai.1$V7 > fai.2$V7[i]] + 1
# 		} else {
# 			system(paste0('sed -i "',fai.2$V7[i],'s/$/N/" fasta/all_species/cds_corrected_all_codons.fna'))
# 		}
		system(paste0('sed -i "',fai.2$V7[i],'s/$/N/" fasta/all_species/cds_corrected_all_codons.fna'))
	}
}

# If the final codon has one nucleotide, pal2nal seems to be okay with it.
# For completeness's sake, remove the last base
if (nrow(fai.1)) {
	for (i in 1:nrow(fai.1)) {
		message('Fixing CDS sequence ',fai.1$V1[i],' by removing one final base.')
		system(paste0('sed -i "',fai.1$V7[i],'s/.$//" fasta/all_species/cds_corrected_all_codons.fna'))
	}
}