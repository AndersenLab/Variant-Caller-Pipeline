library(reshape)
library(compare)
library(plyr)

setwd("../raw/mmp")

load_vcf <- function(infile) {
  mmp <- read.table(infile, header=T,sep="\t", comment.char = "")
  # Fix Names
  names(mmp)[1:2] <- c("chr", "pos")
  # Remove REF Column
  mmp <- subset(mmp, select = -c(REF))
  r <- melt(mmp, id.vars=c("chr","pos"))
  names(r)[3:4] <- c("strain","a_mmp")
  r
}

mmp_snps <- read.table("mmp_snps.txt", header=T, sep="\t", comment.char = "")
names(mmp_snps)[c(1,4)] <- c("chr","t_mmp")


a_mmp <- load_vcf("andersen_mmp_Q10.txt")


m <- join(x=a_mmp,y=mmp_snps, type="left", by = c("chr","pos","strain"))
  
