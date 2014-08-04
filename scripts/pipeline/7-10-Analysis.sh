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
set Qfiles=`ls -l ../../data/vcf/12_new_caller.filter.Q*.bcf | sed 's#^.*/##'`

sbatch 03_vcf_compare.R "Quality Comparison - RadSeq" andersen08_radseq.ws220.bcf 12_new_caller.filter.lcr.bcf $Qfiles
sbatch 03_vcf_compare.R "Quality Comparison - MMP" mmp.vcf.gz 12_new_caller.filter.lcr.bcf $Qfiles
sbatch 03_vcf_compare.R "Depth Comparison - Rockman"  Rockman_ws220.bcf $Qfiles

# Compare Across Different depths

foreach i (100 250 300 350 400 500 avg avg3)
	sbatch --nodelist=node2,node3 02_filters.py 12_new_caller.bcf  -l -d  $i
end

set Dfiles=`ls -l ../../data/vcf/12_new_caller.filter.d*.bcf | sed 's#^.*/##'`

sbatch 03_vcf_compare.R "Depth Comparison - RadSeq" andersen08_radseq.ws220.bcf 12_new_caller.filter.lcr.bcf $Dfiles
sbatch 03_vcf_compare.R "Depth Comparison - MMP" mmp.vcf.gz 12_new_caller.filter.lcr.bcf $Dfiles
sbatch 03_vcf_compare.R "Depth Comparison - Rockman"  Rockman_ws220.bcf $Dfiles