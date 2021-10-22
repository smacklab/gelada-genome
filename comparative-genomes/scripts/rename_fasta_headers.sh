#!/bin/sh

i=$(sed -n ${1}p id_lists/orthogroups/sco_all.txt)

k=0
for j in $(seq $(wc -l assembly_list.txt | cut -d ' ' -f 1)); do
proteins=$(cat id_lists/orthogroups/orthogroups.tsv | grep ${i} | cut -f $(($j + 1)))
if [ $(echo $proteins | grep -c '^$') -gt 0 ] || [ $(echo $proteins | grep -c ',') -gt 0 ]; then
name=$(sort -k1 assembly_list.txt | sed -n ${j}p | cut -f 4)
echo Orthogroup ${i} does not contain an ortholog for taxon ${name}
else
k=$(($k+1))
name=$(sort -k1 assembly_list.txt | sed -n ${j}p | cut -f 4)
sed -i 's/>'${k}'$/>'${name}'/g' fasta/sco/aa/${i}.faa
sed -i 's/>'${k}'$/>'${name}'/g' fasta/sco/cds/${i}.fna
sed -i 's/>'${k}'$/>'${name}'/g' fasta/sco/cds_corrected/${i}.fna
fi
done

exit
