
module load iqtree/2.1.2

# iqtree2 -B 10000 -s cytb_all.aln.phy -o NC_005943 -p data/cytb_partitions.txt
iqtree2  --redo -B 10000 -s cytb_haplotypes.aln.phy -o NC_005943 -m MFP -p data/cytb_partitions.txt
# iqtree2 -B 10000 -s cytb_haplotypes.aln.phy -o FIL001 -p data/cytb_partitions.txt
