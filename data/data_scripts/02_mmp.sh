#!/bin/bash

cd ../raw/mmp

# Download Data from "The million mutation project: A new approach to genetics in Caenorhabditis elegans"
wget 'http://genome.sfu.ca/mmp/mmp_wild_isolate_data_Mar13.txt'

# Generate gff files for each strain
cut -f 3 mmp_wild_isolate_data_Mar13.txt | sort | uniq > strains.txt

egrep -v '(insertion|deletion|complex_change)' mmp_wild_isolate_data_Mar13.txt > mmp_snps.txt

# Generate individual files for each strain (will be concatenated later)
mkdir strains

# Generate unique positions
cut -f 4,5,6 mmp_snps.txt | sort | uniq | sed -e 's/^/chr/' > chr_pos.txt

cat strains.txt | xargs -I {} -n 1 -P 5 sh -c "grep '{}' mmp_snps.txt | cut -f 4,5,6,7 | sed -e 's/^/chr/'  | sort -n -k 1,2 > strains/{}.txt"

join -1 1 -2 2 AB1.txt AB3.txt

# Get variants
gunzip  -kfc 04_mmp_strains.txt.Q10.vcf.gz | vcf-to-tab > mmpQ10.txt

# Convert file to a 'gff'
echo "##fileformat=VCFv4.1" > mmp235.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t`cat strains.txt`" >> mmp235.vcf                      

awk -F $'\t' '{ print $4 "\t" $5 "\t.\T" $6 "\t" $7 "\t.\t.\t.\t" }'  <(tail +2 mmp_wild_isolate_data_Mar13.txt | head -n 50) >> mmp235.vcf
