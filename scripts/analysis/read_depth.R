# Read Depth Stats

bam_depth <- read.table("~/Documents/Spring Rotation/Variant-Caller-Pipeline/data/ancillary/bam_depth.txt", header=T, quote="\"")
View(bam_depth)

vars <- c("Run","Library","Strain","junk")
bam_depth[,vars] <- colsplit(bam_depth$bam, "-|-", vars)
bam_depth$junk <- NULL

fastq_info <- read.csv("~/Documents/Spring Rotation/Variant-Caller-Pipeline/data/ancillary/fastq_info.csv")