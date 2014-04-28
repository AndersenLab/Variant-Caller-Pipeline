#!/bin/bash

# Construct variant calling comparisons

#####
# 1 # ALL Together
#####


sbatch 01_mpileup.py "00_all_bams.txt" 


sbatch 01_mpileup.py "00_all_bams.txt"
sbatch 01_mpileup.py "01a_BGI1_set.txt"
sbatch 01_mpileup.py "01b_BGI2_set.txt"
sbatch 01_mpileup.py "01c_BGI3_set.txt"
sbatch 01_mpileup.py "02a_BGI2_rep1.txt"
sbatch 01_mpileup.py "02b_BGI2_rep2.txt"
sbatch 01_mpileup.py "03_RET1a.txt"
sbatch 01_mpileup.py "03_RET1b.txt"
sbatch 01_mpileup.py "03_RET1.txt"
sbatch 01_mpileup.py "03_RET2a.txt"
sbatch 01_mpileup.py "03_RET2b.txt"
sbatch 01_mpileup.py "03_RET2.txt"
sbatch 01_mpileup.py "03_RET3a.txt"
sbatch 01_mpileup.py "03_RET3b.txt"
sbatch 01_mpileup.py "03_RET3.txt"
sbatch 01_mpileup.py "03_RET4a.txt"
sbatch 01_mpileup.py "03_RET4b.txt"
sbatch 01_mpileup.py "03_RET4.txt"
sbatch 01_mpileup.py "03_RET5a.txt"
sbatch 01_mpileup.py "03_RET5b.txt"
sbatch 01_mpileup.py "03_RET5.txt"
sbatch 01_mpileup.py "03_RET6a.txt"
sbatch 01_mpileup.py "03_RET6b.txt"
sbatch 01_mpileup.py "03_RET6.txt"
sbatch 01_mpileup.py "03_RET7a.txt"
sbatch 01_mpileup.py "03_RET7b.txt"
sbatch 01_mpileup.py "03_RET7.txt"
sbatch 01_mpileup.py "04_mmp_strains.txt"

#####
# 2 # Compare runs indiv. vs. all together
#####

LIBS=`cat "../comparisons/02a_btwn_all_[all].txt" | egrep -o 'RET[0-9][ab]?' | sort | uniq`

for l in $LIBS;
do
	cat "../comparisons/02a_btw_all_[all].txt" | egrep "$l" > ../comparisons/02a_btwn_all_[$l].txt
done

for fset in `ls ../comparisons/02a_btwn_all*`;
do
	sbatch 01_mpileup.py $fset
done

#####
# 3 # BGI2 Repeats
#####

sbatch 01_mpileup.py "../comparisons/02a_BGI2_rep1.txt"
sbatch 01_mpileup.py "../comparisons/03a_BGI2_repeats_[set2].txt"

###################
# Run Comparisons #
###################

sbatch 02_plot_stats.py mem_vs_aln ../comparisons/01a_mem_vs_aln_[mem].txt.vcf ../comparisons/01b_mem_vs_aln_[aln].txt.vcf
plot-vcfstats ../reports/mem_vs_aln.txt -p mem_vs_aln/