# Generate Depth Statistics for a given VCF or set of VCFs

library(ggplot2)
library(reshape2)

#setwd("/Users/dancook/Documents/Spring_Rotation/Variant-Caller-Pipeline/data/vcf/")

#vcfs <-  sprintf("%s", commandArgs(trailingOnly = TRUE))

vcfs = c("00_all_bams.txt.vcf.gz", "00_all_bams.txt.Q40.vcf.gz", "01a_BGI1_set.txt.Q40.vcf.gz", "01b_BGI2_set.txt.Q40.vcf.gz", "01c_BGI3_set.txt.Q40.vcf.gz")

#vcfs = c("03_RET2a.txt.Q40.vcf.gz", "03_RET2b.txt.Q40.vcf.gz")


d <- lapply(vcfs, function(v) { cbind("vcf"=v, "depth"= as.numeric(system(sprintf("gunzip -kfc %s | egrep -o 'DP=[0-9]*' | sed -e 's/DP=//g' ", v), intern=T))) })
names(d) <- vcfs
d <- as.data.frame(do.call("rbind",d))
d$depth <- as.numeric(d$depth)

d <- group_by(d, vcf) %.%
      mutate(mean_depth=mean(depth)) %.%
      mutate(d_plus_2sqrt_d=mean_depth + (2*sqrt(mean_depth))) %.%
      mutate(d_plus_3sqrt_d=mean_depth + (3*sqrt(mean_depth)))

ggplot(d, aes(x=depth, color = vcf)) + 
  stat_ecdf() + 
  labs(title="Cumulative Distribution Frequences" ) +
  geom_vline(aes(xintercept=d_plus_2sqrt_d, color=vcf, linetype = "longdash")) +
  geom_vline(aes(xintercept=d_plus_3sqrt_d, color=vcf, linetype = "dotted"))
  


