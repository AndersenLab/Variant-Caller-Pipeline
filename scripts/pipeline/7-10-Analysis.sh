# Call variants using old/new caller (w/ samtools 2-rc9) and old version of samtools 0.1.19-44428cd
prepare samtools

#set dep=`squeue | grep 'vcfFilt' | awk 'ORS=""; {printf(":%s", $1)}'`

sbatch 03_vcf_compare.R "Caller Comparison" 12_old_samtools.txt.vcf.gz 12_old_caller.bcf 12_new_caller.bcf
sbatch 03_vcf_compare.R "Caller Comparison - No LCR" 12_old_samtools.txt.filter.lcr.bcf 12_old_caller.filter.lcr.bcf 12_new_caller.filter.lcr.bcf 
sbatch 03_vcf_compare.R "Caller Comparison - No LCR and multiallelic removed" 12_old_samtools.txt.filter.lcr.bcf 12_old_caller.filter.lcr.bcf 12_new_caller.filter.lcr.m.bcf 

# Compare across different qualities and depths
foreach i (20 30 40 50 75 100 125 150 175 200 220 240 260 280 300 400 450 500)
	sbatch --nodelist=node2,node3 02_filters.py 12_new_caller.bcf  -l -Q $i
end
#sbatch --dependency="afterok$dep" 03_vcf_compare.R  12_old_caller.filter.lcr.bcf 12_new_caller.filter.lcr.bcf 12_old_samtools.txt.filter.lcr.bcf
set Qfiles=`ls ../../data/vcf/12_new_caller.filter.Q*.bcf | sed 's#^.*/##'`

sbatch 03_vcf_compare.R "Quality Comparison - RadSeq" andersen08_radseq.ws225.bcf 12_new_caller.filter.lcr.bcf $Qfiles
sbatch 03_vcf_compare.R "Quality Comparison - MMP" mmp.vcf.gz 12_new_caller.filter.lcr.bcf $Qfiles


# Compare Across Different depths

sbatch 02_filters.py 12_new_caller.bcf  -l -d 100
sbatch 02_filters.py 12_new_caller.bcf  -l -d 250
sbatch 02_filters.py 12_new_caller.bcf  -l -d 300
sbatch 02_filters.py 12_new_caller.bcf  -l -d 350
sbatch 02_filters.py 12_new_caller.bcf  -l -d 400
sbatch 02_filters.py 12_new_caller.bcf  -l -d 500
sbatch 02_filters.py 12_new_caller.bcf  -l -d 750
sbatch 02_filters.py 12_new_caller.bcf  -l -d avg2
sbatch 02_filters.py 12_new_caller.bcf  -l -d avg3


set Dfiles=`ls ../../data/vcf/12_new_caller.filter.d*.bcf | sed 's#^.*/##'`

sbatch 03_vcf_compare.R "Depth Comparison - RadSeq" andersen08_radseq.ws225.bcf 12_new_caller.filter.lcr.bcf $Dfiles
sbatch 03_vcf_compare.R "Depth Comparison - MMP" mmp.vcf.gz 12_new_caller.filter.lcr.bcf $Dfiles