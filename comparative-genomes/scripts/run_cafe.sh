#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-type=END
#SBATCH --mail-user=kchiou@uw.edu
#SBATCH --job-name="cafe"
#SBATCH --account=csde
#SBATCH --partition=csde
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --mem=8GB

mkdir -p results
mkdir -p results/cafe

mkdir -p id_lists/taxon_joins

cat data/cafe.sh | sed s/chronogram/$(cat data/tree.nwk | sed s/\;//g)/g > results/cafe.sh
chmod +x results/cafe.sh

# /gscratch/scrubbed/smacklab/comparative_genomes/id_lists/orthogroups/orthogroups.tsv

rm -f id_lists/taxon_joins/prot2taxon.txt

for i in $(seq $(wc -l assembly_list.txt | cut -d ' ' -f 1)); do
	acc=$(sed -n ${i}p assembly_list.txt | cut -f 1 | sed "s/ /_/g")
	nam=$(sed -n ${i}p assembly_list.txt | cut -f 2 | sed "s/ /_/g")
	tax=$(sed -n ${i}p assembly_list.txt | cut -f 4)
	protlist=id_lists/feature_joins/${acc}_${nam}_prot2gene.txt
	protcount=$(cat $protlist | wc -l)
	paste <(printf $tax$'\n'%.0s $(seq ${protcount})) <(cut -f 2 id_lists/feature_joins/${acc}_${nam}_prot2gene.txt) >> id_lists/taxon_joins/prot2taxon.txt
done

module load contrib/r/3.6.1

scripts/calculate_orthogroup_sizes.R

module load contrib/cafe/4.2

# caferror.py -d ./ -l logs/caferror.log -o caferror_output.txt -i results/cafe.sh
# mv cafe_* results/cafe
# mv caferror.sh results

caferror.py -d cafe-tmp -i results/cafe.sh

# Rerun CAFE to export p-values for all families, (default is to hide families with non-significant family-wise p-values)
cat caferror.sh | sed 's/ -sp .*//g' | uniq | sed 's/\(^errormodel .*\)/\1 -all/g' | sed 's/\(^load .*\)/\1 -p 1/g' | sed 's/cafe_final/cafe_error_final/g' | sed 's/^lambda -s/lambda -l '$(grep '^Lambda' cafe-tmp/cafe_final_log.txt | tail -n 1 | cut -d ' ' -f 3)'/g' > caferror_final.sh
chmod +x caferror_final.sh

./caferror_final.sh

cp cafe-tmp/cafe_error_final_report.cafe results/cafe
cp caferror_final.sh results/cafe

# module load contrib/cafe/5.0
# cafexp -t data/tree.nwk -i id_lists/orthogroups/protein_counts.txt -p -e $(grep '^errormodel ' caferror_final.sh | sed 's/^errormodel -model \([a-z0-9/._-]*\) -all/\1/')

# cafexp -t data/tree.nwk -i id_lists/orthogroups/protein_counts.txt -p -e results/cafe/cafe_errormodel_0.03173828125.txt

# results/cafe.sh

exit
