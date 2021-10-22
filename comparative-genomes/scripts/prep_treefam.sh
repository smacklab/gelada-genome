#!/bin/bash

git clone https://github.com/treefam/treefam_tools.git
cd treefam_tools/treefam_scan/

mkdir hmm_lib
cd hmm_lib/
wget http://www.treefam.org/static/download/treefam9.hmm3.tar.gz
tar -xzvf treefam9.hmm3.tar.gz

module load contrib/perl/5.26.2
module load contrib/hmmer/3.0

mv treefam9.hmm3 TreeFam9
hmmpress TreeFam9

cd ..

mkdir -p fasta
cp ../../fasta/aa_longest_isoform/*.faa fasta
