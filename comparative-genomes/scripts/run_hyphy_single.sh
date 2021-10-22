#!/bin/bash

wd=$(pwd)

i=$(( ( $1 + 1 ) / 2 ))

orthogroup=$(sed -n ${i}p id_lists/orthogroups/sco_most.txt)

mkdir -p nexus
mkdir -p nexus/cds

if [ ! -f nexus/cds/${orthogroup}.cds.nex ]; then
#	module load contrib/prank/170427
	module load gcc/9.2.0
	export PATH=$PATH:/scratch/klchiou/sw/prank/170427/bin
	prank -convert -d=paml/${orthogroup}/${orthogroup}.paml -f=nexus -o=nexus/cds/${orthogroup}.cds
fi

if [ ! -f hyphy/trees/${orthogroup}.nwk ]; then
#	module load contrib/r/3.6.1
	module load r/latest
	scripts/create_hyphy_tree.R $orthogroup
#	module unload contrib/r/3.6.1
	module unload r/latest
fi

# module load contrib/hyphy/2.5.1
module load hyphy/2.5.9

msa=nexus/cds/${orthogroup}.cds.nex
nwk=hyphy/trees/${orthogroup}.nwk

mkdir -p hyphy
mkdir -p hyphy/busted
mkdir -p hyphy/relax

if [ $(( $1 % 2 )) -eq 1 ]; then
	if [ ! -f hyphy/busted/${orthogroup}.busted.txt ]; then
		hyphy busted --alignment $msa --tree $nwk --branches Foreground --output hyphy/busted/${orthogroup}.json > hyphy/busted/${orthogroup}.out
		#paste <(echo $orthogroup) <(grep '^Likelihood ratio test' hyphy/busted/${orthogroup}.out | tail -n 1 | cut -d '=' -f 2 | sed 's/ *\([0-9.e-]*\).*/\1/g') > hyphy/busted/${orthogroup}.busted.txt
		paste <(echo $orthogroup) <(cat hyphy/busted/${orthogroup}.json | grep '"LRT"' | sed 's/^ *"LRT"://g' | sed 's/,$//g') <(cat hyphy/busted/${orthogroup}.json | grep '"p-value"' | sed 's/^ *"p-value"://g') > hyphy/busted/${orthogroup}.busted.txt
	fi
else
	if [ ! -f hyphy/relax/${orthogroup}.relax.txt ]; then
		hyphy relax --alignment $msa --tree $nwk --test Foreground --output hyphy/relax/${orthogroup}.json > hyphy/relax/${orthogroup}.out
		#paste <(echo $orthogroup) <(grep '^Likelihood ratio test' hyphy/relax/${orthogroup}.out | tail -n 1 | cut -d '=' -f 2 | sed 's/ *\([0-9.e-]*\).*/\1/g') > hyphy/relax/${orthogroup}.relax.txt
		paste <(echo $orthogroup) <(cat hyphy/relax/${orthogroup}.json | grep '"LRT"' | sed 's/^ *"LRT"://g' | sed 's/,$//g') <(cat hyphy/relax/${orthogroup}.json | grep '"p-value"' | sed 's/^ *"p-value"://g') > hyphy/busted/${orthogroup}.relax.txt
	fi
fi

exit