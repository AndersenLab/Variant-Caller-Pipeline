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
  arrange(ext) %.%
  mutate(bam=sprintf("%s-%s-%s-%s.bam", run, library, strain, paste0(unique(hash), collapse='-'))) %.%
  filter(pair==1) %.%
  ungroup() %.%
  select(bam, lane, flowcell)


d <- inner_join(fastq_info, bam_depth, by = 'bam')
d <- rename(d,c("Library"="lib"))

d$lib <- as.factor(d$lib)
d$run <- as.factor(d$Run)

# Plot avg. depth by Run
ggplot(data=d) +
  geom_boxplot( aes(y=avg_depth,x=factor(lib), fill=Run), drop=F) +
  scale_fill_brewer("RdGy",palette="Set4")

d <- group_by(d, Strain) %.%
     mutate(strain_mean=sum(avg_depth)) %.%
     arrange(strain_mean)
# Plot avg. depth by Strain
ggplot(data=d) +
  geom_line( aes(y=strain_mean,x=reorder(Strain,strain_mean), group=1)) +
  scale_y_continuous(breaks=seq(0, 100, 10)) +
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
  labs(x="Strain", y="Depth (Run 1 + Run 2 + Run 3)", title="Average Sequencing Depth by Strain") 

ggsave("avg_depth_by_strain.pdf")
