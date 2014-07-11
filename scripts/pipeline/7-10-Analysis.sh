

# Run A number of Analysis

# Call variants using old/new caller (w/ samtools 2-rc9) and old version of samtools 0.1.19-44428cd
prepare samtools
sbatch --nodelist=node2 01_mpileup.py  12_new_caller.txt
sbatch --nodelist=node2 01_mpileup_c.py 12_old_caller.txt

prepare samtools-0.1.19
prepare vcftools
sbatch --nodelist=node4 01_mpileup_0.1.19.py 12_old_samtools.txt


# Compare with LCR Regions Removed
sbatch 02_filters.py 12_old_caller.bcf -l
sbatch 02_filters.py 12_new_caller.bcf  -l
sbatch 02_filters.py 12_old_samtools.txt.vcf.gz -l

# Remove multiallelic sites
sbatch 02_filters.py 12_new_caller.bcf  -l -m

# Remove multiallelic variants

#set dep=`squeue | grep 'vcfFilt' | awk 'ORS=""; {printf(":%s", $1)}'`
#sbatch --dependency="afterok$dep" 03_vcf_compare.R  12_old_caller.filter.lcr.bcf 12_new_caller.filter.lcr.bcf 12_old_samtools.txt.filter.lcr.bcf
sbatch 03_vcf_compare.R "Caller Comparison" 12_old_samtools.txt.vcf.gz 12_old_caller.bcf 12_new_caller.bcf
sbatch 03_vcf_compare.R "Caller Comparison - No LCR" 12_old_samtools.txt.filter.lcr.bcf 12_old_caller.filter.lcr.bcf 12_new_caller.filter.lcr.bcf 
sbatch 03_vcf_compare.R "Caller Comparison - No LCR and multiallelic removed" 12_old_samtools.txt.filter.lcr.bcf 12_old_caller.filter.lcr.bcf 12_new_caller.filter.lcr.m.bcf 


# Compare across different qualities and depths
sbatch 02_filters.py 12_new_caller.bcf  -l -Q 20
sbatch 02_filters.py 12_new_caller.bcf  -l -Q 30
sbatch 02_filters.py 12_new_caller.bcf  -l -Q 40
sbatch 02_filters.py 12_new_caller.bcf  -l -Q 75
sbatch 02_filters.py 12_new_caller.bcf  -l -Q 100
sbatch 02_filters.py 12_new_caller.bcf  -l -Q 150
sbatch 02_filters.py 12_new_caller.bcf  -l -Q 200
sbatch 02_filters.py 12_new_caller.bcf  -l -Q 300
sbatch 02_filters.py 12_new_caller.bcf  -l -Q 400
sbatch 02_filters.py 12_new_caller.bcf  -l -Q 500

sbatch 03_vcf_compare.R "Quality Comparison - RadSeq" andersen08_radseq.ws225.bcf 12_new_caller.filter.lcr.bcf 12_new_caller.filter.Q0020.lcr.bcf 12_new_caller.filter.Q0030.lcr.bcf 12_new_caller.filter.Q0040.lcr.bcf 12_new_caller.filter.Q0075.lcr.bcf 12_new_caller.filter.Q0100.lcr.bcf 12_new_caller.filter.Q0150.lcr.bcf 12_new_caller.filter.Q0200.lcr.bcf 12_new_caller.filter.Q0300.lcr.bcf 12_new_caller.filter.Q0400.lcr.bcf 12_new_caller.filter.Q0500.lcr.bcf 

sbatch 03_vcf_compare.R "Quality Comparison - MMP" mmp.vcf.gz 12_new_caller.filter.lcr.bcf 12_new_caller.filter.Q0010.lcr.bcf 12_new_caller.filter.Q0020.lcr.bcf 12_new_caller.filter.Q0030.lcr.bcf 12_new_caller.filter.Q0040.lcr.bcf 12_new_caller.filter.Q0050.lcr.bcf 12_new_caller.filter.Q0100.lcr.bcf 


# Compare Across Different depths
sbatch 02_filters.py 12_new_caller.bcf  -l -d 500
sbatch 02_filters.py 12_new_caller.bcf  -l -d avg2
sbatch 02_filters.py 12_new_caller.bcf  -l -d avg3
sbatch 02_filters.py 12_new_caller.bcf  -l -d 1000
sbatch 02_filters.py 12_new_caller.bcf  -l -d 1500
sbatch 02_filters.py 12_new_caller.bcf  -l -d 2000