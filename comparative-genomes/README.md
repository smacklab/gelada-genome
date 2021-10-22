# Gelada comparative genomics pipeline

Scripts in this folder cover analyses performed on the gelada genome in comparison with other genomes, as well as comparisons of gelada morphometric and hematological measurements.

## Comparative genomics analysis

```
scripts/download_genomes.sh

for i in $(seq $(wc -l assembly_list.txt | cut -d ' ' -f 1)); do
scripts/subset_proteins_longest_isoform.sh
done

sbatch scripts/run_orthofinder.sh

scripts/prep_palo.sh
sbatch --array=1-$(wc -l id_lists/orthogroups/sco_all.txt | cut -d ' ' -f 1) scripts/run_palo.sh

scripts/fix_cds_orthogroups.sh

sbatch scripts/rename_fasta_headers_parallel.sh

sbatch --array=1-$(($(wc -l id_lists/orthogroups/sco.txt | cut -d ' ' -f 1) / 100 + 1)) scripts/run_paml_parallel.sh
sbatch --array=1-$(($(wc -l id_lists/orthogroups/sco.txt | cut -d ' ' -f 1) / 50 + 1)) scripts/run_hyphy_parallel.sh

mkdir -p results/paml
cat paml/*/*_out.paml > results/paml/paml_results.txt


# Calculate orthogroup alignment lengths
for i in `cat id_lists/orthogroups/sco_all.txt`; do
echo $i;
paste <(echo $i) <(head -1 paml/${i}/${i}.paml | sed 's/^ *\([0-9]*\) *\([0-9]*\)/\2\t\1/g') >> results/orthogroups_dims_all.txt;
done


rm -f results/orthogroups_dims_mox.txt
for i in `cat id_lists/orthogroups/sco_all.txt`; do
if [ -f paml/${i}/${i}.paml ]; then
echo $i;
paste <(echo $i) <(head -1 paml/${i}/${i}.paml | sed 's/^ *\([0-9]*\) *\([0-9]*\)/\2\t\1/g') >> results/orthogroups_dims_mox.txt;
else
echo $i skipping;
fi
done

rm -f results/orthogroups_dims_agave.txt
for i in `cat id_lists/orthogroups/sco_all.txt`; do
if [ -f paml/${i}/${i}.paml ]; then
echo $i;
paste <(echo $i) <(head -1 paml/${i}/${i}.paml | sed 's/^ *\([0-9]*\) *\([0-9]*\)/\2\t\1/g') >> results/orthogroups_dims_agave.txt;
else
echo $i skipping;
fi
done


# Do treefams
scripts/prep_treefam.sh
scripts/run_treefam.sh
scripts/summarize_treefam.sh

scripts/run_cafe.sh

scripts/summarize_cafe.R
scripts/summarize_paml.R

```

## Morphometric analysis pipeline

```
scripts/lung_volume.R
```

## Hematological analysis pipeline

```
# Hb-O2 analysis
scripts/visualize_p50.R

# Hb concentration analysis
scripts/hb_concentration.R
```