# Read Depth Stats

bam_depth <- read.table("../../data/ancillary/bam_depth.txt", header=T, quote="\"",  strip.white = T)
View(bam_depth)

vars <- c("Run","Library","Strain","junk")
bam_depth[,vars] <- colsplit(bam_depth$bam, "-|-", vars)
bam_depth$junk <- NULL

# Construct bam filename
fastq_info <- read.csv("../../data/ancillary/fastq_info.txt",  strip.white = T)
vars <- c("Run","Library","Strain","hash", "ext")
fastq_info[,vars] <- colsplit(fastq_info$fastq, "-|-", vars)

fastq_info <- group_by(fastq_info,run, library, strain, flowcell, lane) %.%
  mutate(bam=sprintf("%s-%s-%s-%s.bam", run, library, strain, paste0(unique(hash), collapse='-'))) %.%
  filter(pair==1) %.%
  ungroup() %.%
  select(bam, lane, flowcell)

d <- inner_join(fastq_info, bam_depth, by = 'bam')
d <- rename(d,c("Library"="lib"))

ggplot(data=d, aes(y=avg_depth,x=lib, order=-desc(lib))) +
  geom_boxplot() +
  coord_flip() +
  facet_grid(Run ~ .) 


# Format Fastq info
