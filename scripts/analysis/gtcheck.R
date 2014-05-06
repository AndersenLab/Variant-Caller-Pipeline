library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)


args <- c("01a_mem_vs_aln_[mem].txt.vcf.gz","01b_mem_vs_aln_[aln].txt.vcf.gz","01a_mem_vs_aln_[mem].txt.vcf.gz","03a_BGI2_repeats_[set2].txt.vcf.gz")
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

# Generate concordance numbers
concordance <- function(f1, f2) {
  # Create temporary file
  tmp <- tempfile()
  system(sprintf('bcftools gtcheck -s  %s %s > %s', f1, f2, tmp))
  SM <- import_table("SM", tmp)
  CN <- import_table("CN", tmp)
  SM$comparison <- sprintf('%s__%s', sub(".txt","",sub(".vcf.gz","",f1)), sub(".txt","",sub(".vcf.gz","",f2)))
  
  SM$Average.Discordance.Number.of.sites <- runif(length(SM$Average.Discordance.Number.of.sites))
  
  list(SM=SM,CN=CN)
}

results <- lapply(pairs, function(x) { concordance(x[1],x[2]) })

# Generate Concordance dataframe
c_frame <- do.call(rbind.data.frame,sapply(results, function(x) { x["SM"] }))

# Generate strain average for sorting.
c_frame <- group_by(c_frame, Sample) %.% 
  mutate(sample_avg_depth=mean(Average.Discordance.Number.of.sites)) %.%
  arrange(sample_avg_depth)
# Plot concordance of multiple concordance results
ggplot(c_frame) +
    geom_line(position="identity",aes(y=as.numeric(Average.Discordance.Number.of.sites), x=reorder(c_frame$Sample, c_frame$sample_avg_depth), color=comparison, group=comparison)) +
    labs(title="Average Discordance", x="Sample Name", y="Average Discordance (%)") +
    scale_y_continuous(breaks=seq(0, 10, 0.1))

ggsave(file=sprintf("../reports/quality/%s.png", args[1] ))

