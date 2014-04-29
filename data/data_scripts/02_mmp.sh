#!/bin/bash

cd ../raw/mmp

# Download Data from "The million mutation project: A new approach to genetics in Caenorhabditis elegans"
wget 'http://genome.sfu.ca/mmp/mmp_wild_isolate_data_Mar13.txt'

# Generate individual files for each strain (will be concatenated later)
mkdir strains

# Generate set of strains
cut -f 3 mmp_wild_isolate_data_Mar13.txt | sort | uniq > strains.txt

# Generate list of SNPs , filtering out indels & complex changes.
egrep -v '(insertion|deletion|complex_change)' mmp_wild_isolate_data_Mar13.txt | awk -F $'\t' '{ print "CHROMOSOME_" $4 "\t" $3 "\t" $6 "/" $7 "\t" $5 "\t" $5+1 "\t0\t+\t" 0 "\t"}' |  tail -n +2  > mmp_snps.txt

perl remap_gff_between_releases.pl --release1 235 --release2 220 --out out.test --gff mmp_snps.txt 


# Perform liftover

cd ../vcf/

gunzip -kfc ../vcf/04_mmp_strains.txt.Q10.vcf.gz | vcf-to-tab > andersen_mmp_Q10.txt


##############
# Plot stats #
##############

bcftools stats  -s - 02a_BGI2_rep1.txt.Q40.vcf.gz 02b_BGI2_rep2.txt.Q40.vcf.gz > rep1.txt
plot-vcfstats -p rep1_rep2/ rep1.txt