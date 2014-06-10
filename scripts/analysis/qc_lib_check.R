# Checking Strains

load("/Users/dancook/Documents/Spring Rotation/Variant-Caller-Pipeline/data/reports/X01a_BGI1_set.txt.Q40.vcf.gzBGI2.RET6.rep1.txt.Q40.vcf.gz.Rdata")

f <- group_by(concordance_results[[1]], Query) %.%
  filter(max(Concordance)==Concordance)


sample_lib <- read.csv("~/Documents/Spring Rotation/Variant-Caller-Pipeline/data/ancillary/sample_lib.csv")
sample_lib$Sample <- sample_lib$Strain
excel(inner_join(f, sample_lib, by="Sample"))

load("/Users/dancook/Documents/Spring Rotation/Variant-Caller-Pipeline/data/reports/X01a_BGI1_set.txt.Q40.vcf.gzBGI2.RET7.rep1.txt.Q40.vcf.gz.Rdata")

f <- group_by(concordance_results[[1]], Query) %.%
  filter(Query=='QX1212')


sample_lib <- read.csv("~/Documents/Spring Rotation/Variant-Caller-Pipeline/data/ancillary/sample_lib.csv")
sample_lib$Sample <- sample_lib$Strain
excel(inner_join(f, sample_lib, by="Sample"))