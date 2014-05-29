#!/bin/bash


# Compare Ret Libraries

Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET1.md')" BGI1-RET1.txt.Q40.vcf.gz BGI2-RET1.txt.Q40.vcf.gz BGI3-RET1.txt.Q40.vcf.gz 
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET2.md')" BGI1-RET2.txt.Q40.vcf.gz BGI2-RET2.txt.Q40.vcf.gz BGI3-RET2.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET3.md')" BGI1-RET3.txt.Q40.vcf.gz BGI2-RET3.txt.Q40.vcf.gz BGI3-RET3.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET4.md')" BGI1-RET4.txt.Q40.vcf.gz BGI2-RET4.txt.Q40.vcf.gz BGI3-RET4.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET5.md')" BGI1-RET5.txt.Q40.vcf.gz BGI2-RET5.txt.Q40.vcf.gz BGI3-RET5.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET6.md')" BGI1-RET6.txt.Q40.vcf.gz BGI2-RET6.txt.Q40.vcf.gz BGI3-RET6.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET7.md')" BGI1-RET7.txt.Q40.vcf.gz BGI2-RET7.txt.Q40.vcf.gz BGI3-RET7.txt.Q40.vcf.gz


Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET1.md')" 03_RET1.txt.Q40
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET2.md')" BGI1-RET2.txt.Q40.vcf.gz BGI2-RET2.txt.Q40.vcf.gz BGI3-RET2.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET3.md')" BGI1-RET3.txt.Q40.vcf.gz BGI2-RET3.txt.Q40.vcf.gz BGI3-RET3.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET4.md')" BGI1-RET4.txt.Q40.vcf.gz BGI2-RET4.txt.Q40.vcf.gz BGI3-RET4.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET5.md')" BGI1-RET5.txt.Q40.vcf.gz BGI2-RET5.txt.Q40.vcf.gz BGI3-RET5.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET6.md')" BGI1-RET6.txt.Q40.vcf.gz BGI2-RET6.txt.Q40.vcf.gz BGI3-RET6.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET7.md')" BGI1-RET7.txt.Q40.vcf.gz BGI2-RET7.txt.Q40.vcf.gz BGI3-RET7.txt.Q40.vcf.gz



# Compare BGI Sets against one another
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI_all.md')" 01a_BGI1_set.txt.Q40.vcf.gz 01b_BGI2_set.txt.Q40.vcf.gz 03c_BGI3_set.txt.Q40.vcf.gz 


# Cleanup Reports folder

gimli # Generates all the pdfs!

rm *.md