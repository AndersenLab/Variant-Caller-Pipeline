
setwd("/Users/dancook/Documents/Spring Rotation/Variant-Caller-Pipeline/data/data_scripts")

library(VariantAnnotation)

vcf <- readVcf("../raw/vcf/test.txt.Q10.vcf.gz", genome="ce10")


# Plots

hist(qual(vcf))