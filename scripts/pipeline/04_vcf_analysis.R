library(reshape2)
library(dplyr)
library(ggplot2)

# Read in vcf file
df <- tbl_df(read.csv2(pipe("bcftools filter -R ../../data/ancillary/ce10_no_LCR.bed test5.bcf | bcftools filter --soft-filter 'QUAL' --include '%QUAL > 30' | egrep -v '(ID=VDB|ID=QS)' | vcf_melt"), as.is=T, sep='\t',header=T, check.names=FALSE))

# Split Prior Likelihoods (PL) of Hom_Ref, Het, Hom_Alt
df <- cbind(df, colsplit(df$PL, ',', names =  c('PL.REF',"PL.HET","PL.ALT")))
df$PL <- NULL

# DP4 Number of alleles used in variant calling:
#    1) forward ref alleles; 
#    2) reverse ref 
#    3) forward non-ref
#    4) reverse non-ref alleles
#    Sum can be smaller than DP because low-quality bases are not counted.

df <- cbind(df, colsplit(df$DP4, ',', names = c('DP4.FORWARD.REF',"DP4.REVERSE.REF","DP4.FORWARD.ALT", "DP4.REVERSE.ALT")))
df$DP4 <- NULL


df$DP4_REF <- df$DP4.FORWARD.REF + df$DP4.REVERSE.REF
df$DP4_ALT <- df$DP4.FORWARD.ALT + df$DP4.REVERSE.ALT


# Split GP Reads = Phred scaled genotype posterior probabilities. Intended for imputed genotype probs. 
df <- cbind(df, colsplit(df$GP, ',', names = c('GP.REF',"GP.HET","GP.ALT")))
df$GP <- NULL

df <- filter(df, DP > 3, DP < 250)

df$GT[df$FILTER == "['QUAL']"] <- "Q"

p <- ggplot(df, aes(DP4_REF, DP4_ALT, color=GT))
p + geom_point(alpha=0.25)


# Aggregate correlation

