

# Compare vcf files!

sbatch 02_filters.py test_old_samtools.txt.vcf.gz -l
sbatch 02_filters.py test6.bcf -l -m
sbatch 02_filters.py test7.bcf -l
sbatch 02_filters.py test2.bcf -l -d avg2 -m
sbatch 02_filters.py test2.bcf -l -d avg3
sbatch 02_filters.py test5.bcf -l -Q 30
sbatch 02_filters.py test2.bcf -l -Q 30 -d 100
set dep=`squeue | grep 'vcfFilt' | awk 'ORS=""; {printf(":%s", $1)}'`

sbatch --dependency="afterok$dep" 03_vcf_compare.R  test_old_samtools.txt.lcr.bcf test_new_caller.lcr.bcf test_old_caller.lcr.bcf
