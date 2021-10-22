#!/bin/bash

cd treefam_tools/treefam_scan

for i in `ls results/*_results.txt`; do
echo 'Now formatting protein families for file '${i}
cat <(echo 'seq_id hmm_name bit_score') <(sed -n 30,$(wc -l $i  | cut -d ' ' -f 1)p $i  | tr -s ' ' | cut -d ' ' -f 1,6,9 | sed 's/[gi|0-9]*ref|\([A-Z0-9_.]*\)| /\1 /g') > results/treefam_hmm_$(echo $i | sed -r 's/results\/treefam_hmm_(.*?)_results.txt/\1/g')_profams.txt
done;

for i in `ls results/*_profams.txt`; do
echo 'Now extracting best protein families for file '${i}
cat <(head -n 1 $i) <(tail -n +2 $i | sort -k1,1 -k3,3rn | sort -uk1,1) > results/treefam_hmm_$(echo $i | sed -r 's/results\/treefam_hmm_(.*?)_profams.txt/\1/g')_bestfams.txt
done;

scripts/summarize_treefam.R

./cafe.R

exit