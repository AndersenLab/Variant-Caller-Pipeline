#!/bin/bash



Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET1.md')" 03_RET1b.txt.Q40.vcf.gz 03_RET1.txt.Q40.vcf.gz 03_RET1a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET2.md')" 03_RET2b.txt.Q40.vcf.gz 03_RET2.txt.Q40.vcf.gz 03_RET2a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET3.md')" 03_RET3b.txt.Q40.vcf.gz 03_RET3.txt.Q40.vcf.gz 03_RET3a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET4.md')" 03_RET4b.txt.Q40.vcf.gz 03_RET4.txt.Q40.vcf.gz 03_RET4a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET5.md')" 03_RET5b.txt.Q40.vcf.gz 03_RET5.txt.Q40.vcf.gz 03_RET5a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET6.md')" 03_RET6b.txt.Q40.vcf.gz 03_RET6.txt.Q40.vcf.gz 03_RET6a.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/RET7.md')" 03_RET7b.txt.Q40.vcf.gz 03_RET7.txt.Q40.vcf.gz 03_RET7a.txt.Q40.vcf.gz

Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1-BGI2-BGI3-RET1.md')" BGI2-RET1.txt.Q40.vcf.gz BGI3-RET1.txt.Q40.vcf.gz 
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1-BGI2-BGI3-RET2.md')" BGI1-RET2.txt.Q40.vcf.gz BGI2-RET2.txt.Q40.vcf.gz BGI3-RET2.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1-BGI2-BGI3-RET3.md')" BGI1-RET3.txt.Q40.vcf.gz BGI2-RET3.txt.Q40.vcf.gz BGI3-RET3.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1-BGI2-BGI3-RET4.md')" BGI1-RET4.txt.Q40.vcf.gz BGI2-RET4.txt.Q40.vcf.gz BGI3-RET4.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1-BGI2-BGI3-RET5.md')" BGI1-RET5.txt.Q40.vcf.gz BGI2-RET5.txt.Q40.vcf.gz BGI3-RET5.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1-BGI2-BGI3-RET6.md')" BGI1-RET6.txt.Q40.vcf.gz BGI2-RET6.txt.Q40.vcf.gz BGI3-RET6.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1-BGI2-BGI3-RET7.md')" BGI1-RET7.txt.Q40.vcf.gz BGI2-RET7.txt.Q40.vcf.gz BGI3-RET7.txt.Q40.vcf.gz



Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI_comp.md')" 01a_BGI1_set.txt.Q40.vcf.gz 01b_BGI2_set.txt.Q40.vcf.gz 03c_BGI3_set.txt.Q40.vcf.gz 01a_BGI1_set.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/all_bgi1.md')" 00_all_bams.txt.Q40.vcf.gz 01a_BGI1_set.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/all_bgi2.md')" 00_all_bams.txt.Q40.vcf.gz 01b_BGI2_set.txt.Q40.vcf.gz 
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/all_bgi3.md')" 00_all_bams.txt.Q40.vcf.gz 01c_BGI3_set.txt.Q40.vcf.gz

Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/rep.md')" 02a_BGI2_rep1.txt.Q40.vcf.gz 02b_BGI2_rep2.txt.Q40.vcf.gz


Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI_test.md')" BGI1-RET1.txt.Q40.vcf.gz  01c_BGI3_set.txt.Q40.vcf.gz

Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI_test2.md')" BGI1-RET1.txt.Q40.vcf.gz  01c_BGI3_set.txt.Q40.vcf.gz




Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1-BGI2-RET1.md')" BGI1-RET1.txt.Q40.vcf.gz BGI2-RET1.txt.Q40.vcf.gz


Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1-BGI2-RET1.md')" BGI1-RET1.txt.Q40.vcf.gz BGI2-RET1.txt.Q40.vcf.gz



Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1-BGI2-RET2.md')" BGI1-RET2.txt.Q40.vcf.gz BGI2-RET2.txt.Q40.vcf.gz
Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1-BGI2-RET3.md')" BGI1-RET3.txt.Q40.vcf.gz BGI2-RET3.txt.Q40.vcf.gz

Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI2_RET6_corrected')" BGI1-RET2.txt.Q40.vcf.gz   BGI2-RET6-rep1.txt.Q40.vcf.gz

Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1_rep1_RET7.md')" 01a_BGI1_set.txt.Q40.vcf.gz  BGI2-RET7-rep1.txt.Q40.vcf.gz

Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI2_RET7_corrected')" BGI1-RET3.txt.Q40.vcf.gz   BGI2-RET7-rep1.txt.Q40.vcf.gz


Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1_rep1_RET4')" 01a_BGI1_set.txt.Q40.vcf.gz   BGI2-RET4-rep1.txt.Q40.vcf.gz



Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../data/reports/BGI1_rep1_RET4')" 01a_BGI1_set.txt.Q40.vcf.gz   BGI2-RET4-rep1.txt.Q40.vcf.gz
