#!/bin/bash



Rscript -e "library(knitr); knit('vcf_report.Rmd', output='../../reports/RET1.md')" 03_RET1b.txt.Q40.vcf.gz 03_RET1.txt.Q40.vcf.gz 03_RET1a.txt.Q40.vcf.gz
