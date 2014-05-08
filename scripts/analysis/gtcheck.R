library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)


args <- c("01a_mem_vs_aln_[mem].txt.vcf.gz","01a_mem_vs_aln_[mem].txt.vcf.gz","03a_BGI2_repeats_[set2].txt.vcf.gz")
pairs <- list()
for (i in 1:(length(args)-1)) {
  pairs<-append(pairs,list(c(args[i], args[i+1])))
}

#---------------------------#
# Import data from bcfstats #
#---------------------------#

# Pipe in gtcheck data from bcfstats
import_table <- function(t_name, f) { 
  # This function can import data from bcftools stats function
  t <- read.csv2(pipe(sprintf("egrep '^(# \\[1\\])?%s[^,]' %s", t_name, f)), as.is=T, sep='\t',header=T, check.names=FALSE)
  names(t) <- lapply(names(t), function(x) make.names(gsub("\\[[1-9]+\\]","",x)))
  t[,c(-1)]
}

ind_results <- function(f) {
  # This function retrieves individual data for each VCF for plotting and comparison.
  tmp <- tempfile()
  system(sprintf('bcftools stats -s -  %s > %s', f, tmp))
  tmp_SN  <- import_table("SN", tmp)
  QUAL <- import_table("QUAL", tmp)
  PSC <- import_table("PSC", tmp)
  
  # Fix summary numbers
  SN <- as.list(tmp_SN[[3]])
  names(SN) <- make.names(tmp_SN[[2]])
  
  list(SN=SN, QUAL=QUAL, PSC=PSC)
}

# Generate concordance numbers
pair_results <- function(f1, f2) {
  # This function retrieves pairwise stats using bcftools for plotting and comparison.
  tmp <- tempfile()
  system(sprintf('bcftools gtcheck -s  %s %s > %s', f1, f2, tmp))
  SM <- import_table("SM", tmp)
  CN <- import_table("CN", tmp)
  SM$comparison <- sprintf('%s__%s', sub(".txt","",sub(".vcf.gz","",f1)), sub(".txt","",sub(".vcf.gz","",f2)))
  
  SM$Average.Discordance.Number.of.sites <- runif(length(SM$Average.Discordance.Number.of.sites))
  
  list(SM=SM,CN=CN)
}

#------------------#
# Individual Stats #
#------------------#

results <- lapply(args, function(x) { ind_results(x) } )

#------------------#
# Pairwise Stats   #
#------------------#

results <- lapply(pairs, function(x) { pair_results(x[1],x[2]) })

#-------------------#
# Concordance Chart #
#-------------------#

# Generate Concordance dataframe
SM_frame <- do.call(rbind.data.frame,sapply(results, function(x) { x["SM"] }))

# Generate strain average for sorting.
SM_frame <- group_by(SM_frame, Sample) %.% 
  mutate(sample_avg_depth=mean(Average.Discordance.Number.of.sites)) %.%
  arrange(sample_avg_depth)

# Plot concordance of multiple concordance results
ggplot(SM_frame) +
    geom_line(position="identity",aes(y=as.numeric(Average.Discordance.Number.of.sites), x=reorder(SM_frame$Sample, SM_frame$sample_avg_depth), color=comparison, group=comparison)) +
    labs(title="Average Discordance", x="Sample Name", y="Average Discordance (%)") +
    scale_y_continuous(breaks=seq(0, 10, 0.1)) +
    theme(legend.position="top")

ggsave(file=sprintf("../reports/quality/%s.png", args[1] ))

#-------------------#
# Concordance Chart #
#-------------------#

con_chart <- function(record) {
  CN_set <- record$CN
  SM_set <- record$SM
  CN_set$Discordance <- runif(length(CN_set$Discordance))
  ggplot(CN_set) +
    geom_tile(aes(x=Sample.i, y=Sample.j, fill=Discordance),colour = "white") + 
    geom_tile(aes(x=Sample.j, y=Sample.i, fill=Discordance),colour = "white") + 
    geom_tile(data=SM_set, aes(x=Sample, y=Sample, fill= Average.Discordance.Number.of.sites)) +
    labs(title=sprintf("%s",SM_set[1,'comparison'])) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme(legend.position="bottom")
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend <- g_legend(con_chart(results[[1]]))

# Generate Concordance Charts and place them within a frame.
png("concordance.png", width = 1200, height=round(length(pairs) / 2)*600, units = "px")
do.call(grid.arrange, c(lapply(results, function(x) { con_chart(x)  }), ncol=2))
dev.off()
