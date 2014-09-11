# Generate Depth Statistics for a given VCF or set of VCFs

library(ggplot2)
library(reshape2)
library(stringr)
# Split out variables

bam_depth <- read.delim("~/Documents/Spring_Rotation/Variant-Caller-Pipeline/data/ancillary/bam_depth.txt")
bam_depth[,c("Run","Library","Strain","hash")] <- str_split_fixed(bam_depth$bam, "-", n=4)

# Box Plot of depths
bam_depth <- group_by(bam_depth, Strain) %>%
  summarise(cum_depth=sum(avg_depth))

# Plot cumulative Depth
ggplot(bam_depth) +
  geom_boxplot(aes(x="Cumulative Depth", y =cum_depth)) +
  theme_bw() +
  labs(x="", y="Whole Genome Sequencing Depth")

ggplot(bam_depth) +
  geom_line(aes(x=Strain, y=cum_depth, group=1), stat="identity", color="#0080ff", size=2) +
  labs(y="Cumulative Depth {avg of several sequencing runs}") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, face="bold", hjust=1))