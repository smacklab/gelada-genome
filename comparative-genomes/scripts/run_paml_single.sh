#!/bin/sh

wd=$(pwd)

# orthogroup=$(sed -n ${1}p id_lists/orthogroups/sco_all.txt)
orthogroup=$(sed -n ${1}p id_lists/orthogroups/sco_most.txt)

if [ ! -f paml/${orthogroup}/${orthogroup}_out.paml ]; then

#	module load contrib/clustalomega/1.2.4
#	module load contrib/perl/5.26.2
#	module load contrib/blast/20021001
#	module load contrib/pal2nal/v14

	module load perl/5.26.0
	export PATH=/scratch/klchiou/sw/clustalomega/1.2.4/bin:/scratch/klchiou/sw/blast/20021001/bin:/scratch/klchiou/sw/pal2nal/14.0:$PATH

	mkdir -p fasta/aln
	mkdir -p fasta/aln/aa

	# Had been removing ambiguous amino acids, but this is actually creating conflicts so reverse course
	# sed -i '/^>/! s/X/-/g' fasta/sco/aa/${orthogroup}.faa
	sed -i '/^>/! s/-/X/g' fasta/sco/aa/${orthogroup}.faa

	clustalo --force -i fasta/sco/aa/${orthogroup}.faa -o fasta/aln/aa/${orthogroup}.aln.faa

	mkdir -p paml
	mkdir -p paml/${orthogroup}

	# Had been removing ambiguous nucleotides, but this is actually creating conflicts so reverse course
	# sed -i '/^>/! s/N/-/g' fasta/sco/cds_corrected/${orthogroup}.fna
	sed -i '/^>/! s/-/N/g' fasta/sco/cds_corrected/${orthogroup}.fna

	# Proteins containing selenocysteine (U) are interpreted as stop codons. Convert these to ambiguous bases.
	sed -i '/^>/! s/U/X/g' fasta/aln/aa/${orthogroup}.aln.faa 

	# Convert protein alignments into codon alignments
	pal2nal.pl fasta/aln/aa/${orthogroup}.aln.faa fasta/sco/cds_corrected/${orthogroup}.fna -output paml -nogap -nomismatch > paml/${orthogroup}/${orthogroup}.paml

	# Set up PAML analysis
	sed -e "s/orthogroup/${orthogroup}/g" data/codeml_BS.H0.ctl > paml/${orthogroup}/${orthogroup}.BS.H0.ctl
	sed -e "s/orthogroup/${orthogroup}/g" data/codeml_BS.H1.ctl > paml/${orthogroup}/${orthogroup}.BS.H1.ctl

#	module load contrib/r/3.6.1
	module load r/latest
	scripts/create_paml_tree.R $orthogroup
#	module load contrib/r/3.6.1
	module load r/latest

#	module load contrib/paml/4.9
	export PATH=/scratch/klchiou/sw/paml/4.9/bin:$PATH

	cd paml/${orthogroup}

	if [ $(grep -c '#1' tree.txt) -gt 0 ]; then

		if [[ -f ${orthogroup}.bs2.h0.txt && $(tail -n 1 ${orthogroup}.bs2.h0.txt | cut -d ':' -f 1) == 'Time used' ]]; then
			echo "Skipping H0 model"
		else
			codeml ${orthogroup}.BS.H0.ctl
		fi

		if [[ -f ${orthogroup}.bs2.h1.txt && $(tail -n 1 ${orthogroup}.bs2.h1.txt | cut -d ':' -f 1) == 'Time used' ]]; then
			echo "Skipping H1 model"
		else
			codeml ${orthogroup}.BS.H1.ctl
		fi

		h0=$(grep lnL ${orthogroup}.bs2.h0.txt | tr -s ' ' | cut -d ' ' -f 5)
		h1=$(grep lnL ${orthogroup}.bs2.h1.txt | tr -s ' ' | cut -d ' ' -f 5)
		T=$(echo "2 * (($h1)-($h0))" | bc)
		p=$(chi2 1 $T | grep prob | tr -s ' ' | cut -d ' ' -f 6)

		echo ${orthogroup}$'\t'${h0}$'\t'${h1}$'\t'${T}$'\t'${p} > ${orthogroup}_out.paml
	
	else

		touch ${orthogroup}_out.paml

	fi

fi

exit
