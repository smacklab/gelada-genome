#!/bin/bash

module load phast/1.3

region=`sed -n ${1}p mammal_sco.txt`

cat mammal_sco/${ortho}.pal2nal.fa | sed 's:[RYSWKMBDHVN]:-:g' > mammal_sco_unambiguous/${ortho}.fa

phyloFit --tree mammal_tree.nwk --out-root=fit/${ortho} mammal_sco_unambiguous/${ortho}.fa

phyloP --subtree=gelada --method=SPH fit/${ortho}.mod mammal_sco_unambiguous/${ortho}.fa > phylop/${ortho}.txt

