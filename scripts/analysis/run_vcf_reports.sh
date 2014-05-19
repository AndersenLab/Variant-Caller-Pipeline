#!/bin/bash



Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET1.md')" 03_RET1b.txt.Q40.vcf.gz 03_RET1.txt.Q40.vcf.gz 03_RET1a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET2.md')" 03_RET2b.txt.Q40.vcf.gz 03_RET2.txt.Q40.vcf.gz 03_RET2a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET3.md')" 03_RET3b.txt.Q40.vcf.gz 03_RET3.txt.Q40.vcf.gz 03_RET3a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET4.md')" 03_RET4b.txt.Q40.vcf.gz 03_RET4.txt.Q40.vcf.gz 03_RET4a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET5.md')" 03_RET5b.txt.Q40.vcf.gz 03_RET5.txt.Q40.vcf.gz 03_RET5a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET6.md')" 03_RET6b.txt.Q40.vcf.gz 03_RET6.txt.Q40.vcf.gz 03_RET6a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET7.md')" 03_RET7b.txt.Q40.vcf.gz 03_RET7.txt.Q40.vcf.gz 03_RET7a.txt.Q40.vcf.gz


Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI_comp.md')" 01a_BGI1_set.txt.Q40.vcf.gz 01b_BGI2_set.txt.Q40.vcf.gz 03c_BGI3_set.txt.Q40.vcf.gz 01a_BGI1_set.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/all_bgi1.md')" 00_all_bams.txt.Q40.vcf.gz 01a_BGI1_set.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/all_bgi2.md')" 00_all_bams.txt.Q40.vcf.gz 01b_BGI2_set.txt.Q40.vcf.gz 
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/all_bgi3.md')" 00_all_bams.txt.Q40.vcf.gz 01c_BGI3_set.txt.Q40.vcf.gz