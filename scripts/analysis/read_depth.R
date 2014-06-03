# Read Depth Stats

bam_depth <- read.table("../../data/ancillary/bam_depth.txt", header=T, quote="\"",  strip.white = T)
View(bam_depth)

vars <- c("Run","Library","Strain","junk")
bam_depth[,vars] <- colsplit(bam_depth$bam, "-|-", vars)
bam_depth$junk <- NULL

bam_depth <- group_by(bam_depth, Strain) %>%
            mutate(cum_depth=sum(avg_depth))

# Cumulative Depth
ggplot(data=bam_depth[bam_depth$Strain!="",]) +
  geom_line( aes(y=cum_depth,x=reorder(Strain,cum_depth), group=1)) +
  scale_y_continuous(breaks=seq(0, 200, 10)) +
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
  labs(x="Strain", y="Depth (Run 1 + Run 2 + Run 3)", title="Cumulative Sequencing Depth by Strain")

ggsave("../../results/cumulative_depth_by_strain.pdf", width=20)


# Depth by Run
ggplot(data=bam_depth[bam_depth$Strain!="",]) +
  geom_line( aes(y=avg_depth,x=reorder(Strain,avg_depth), group=Run, color=Run)) +
  scale_y_continuous(breaks=seq(0, 120, 10), limits=c(0,100)) +
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
  labs(x="Strain", y="Average Depth)", title="Average Sequencing Depth by Strain by Run")

ggsave("../../results/avg_depth_by_Run.pdf", width=20)

# Depth by Library
ggplot(data=bam_depth[bam_depth$Strain!="",]) +
  geom_line( aes(y=avg_depth, x=Strain, group=Library, color=Library), stat="identity") +
  scale_y_continuous(breaks=seq(0, 120, 10), limits=c(0,100)) +
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, size=8)) +
  labs(x="Library", y="Average Depth", title="Average Sequencing Depth by Library") +
  facet_grid(. ~ Library, scales="free_x") +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

ggsave("../../results/avg_depth_by_library.pdf", width=20)
