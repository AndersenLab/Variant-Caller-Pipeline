library(ggplot2)
library(gridExtra)
library(scales)

setwd("~/Documents/Spring Rotation/Variant-Caller-Pipeline/data/vcf/")

args <- commandArgs(trailingOnly = TRUE)

args <- c("03_RET7.txt.Q40.vcf.gz","03_RET7a.txt.Q40.vcf.gz","03_RET7a.txt.Q40.vcf.gz","03_RET7.txt.Q40.vcf.gz")
pairs <- list()
for (i in 1:(length(args)-1)) {
  pairs<-append(pairs,list(c(args[i], args[i+1])))
}

#---------------------------------------#
# Define Functions for generating stats #
#---------------------------------------#

# Function used for plotting later.
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# Pipe in gtcheck data from bcfstats
import_table <- function(t_name, f) { 
  # This function can import data from bcftools stats function
  t <- read.csv2(pipe(sprintf("egrep '^(# )%s[^,]|^%s' %s", t_name, t_name, f)), as.is=T, sep='\t',header=T, check.names=FALSE)
  names(t) <- lapply(names(t), function(x) make.names(gsub("\\[[1-9]+\\]","",x)))
  t[,c(-1)]
}

# Generate a concordance Matrix
get_con_matrix <- function(f1, f2) {
  tmp <- tempfile()
  # Get sample names from each file
  f1_samples <- unlist(str_split(system(sprintf("bcftools view -h %s| grep '#CHROM' | cut -f '10-'",f1), intern=T),'\t'))
  f2_samples <- unlist(str_split(system(sprintf("bcftools view -h %s| grep '#CHROM' | cut -f '10-'",f2), intern=T),'\t'))
  # Get intersection of samples
  joint_samples <- intersect(f1_samples, f2_samples)
  # Retrieve Concordance Data
  r <- read.csv2(pipe(sprintf("echo %s | xargs -I {} -n 1 -P 50 sh -c \"bcftools gtcheck -S {} -g %s %s | sed \'s/$/\t{}/\'\" | egrep -v '#' ", paste0(joint_samples, collapse=" "), f1, f2)), as.is=T, sep='\t',header=F, col.names= c("Concordance", "Uncertainty", "Average Depth", "Number of Sites", "Sample", "Query"))
  r$Comparison <- sprintf('%s__%s.txt', sub(".txt","",sub(".vcf.gz","",f1)), sub(".txt","",sub(".vcf.gz","",f2)))
  r
}

get_pair_results <- function(f1, f2) {
  # This function retrieves individual data for each VCF for plotting and comparison.
  tmp <- 'tmp.txt'
  system(sprintf('bcftools stats -c both -s -  %s %s > %s', f1, f2, tmp))
  SN  <- import_table("SN", tmp)
  QUAL <- import_table("QUAL", tmp)
  PSC <- import_table("PSC", tmp)
  
  comp <- sprintf('%s__%s.txt', sub(".txt","",sub(".vcf.gz","",f1)), sub(".txt","",sub(".vcf.gz","",f2)))
  id_repl <- list('0'=f1, '1'=f2, '2'=comp)
  
  # Replace ID Field as is appropriate
  SN$id <- id_repl[as.character(SN$id)]
  SN$key <- tolower(SN$key)
  QUAL$id <- id_repl[as.character(QUAL$id)]
  PSC$id <- id_repl[as.character(PSC$id)]
  
  # Fix SN Group
  SN <- reshape(SN, direction="wide", timevar="key", ids="id")
  names(SN) <- make.names(gsub("value.","",gsub(":","", names(SN))))
  
  # Add Comparison Column
  SN$Comparison <- comp
  QUAL$Comparison <- comp
  PSC$Comparison <- comp
    
  list(SN=SN, QUAL=QUAL, PSC=PSC)
}


f<-reshape(SN, direction="wide", timevar="key", ids="id")

s <- melt(SN, id.vars="id")

#---------------#
# Generate Data #
#---------------#

pair_results <- lapply(pairs, function(x) { get_pair_results(x[1],x[2]) })

concordance_results <- lapply(pairs, function(x) { get_con_matrix(x[1],x[2]) })


#-----------------------------------------------#
# Plot Individual vcf results for SNPs + Indels #
#-----------------------------------------------#

snp_indel <- cbind(filename=args,melt(do.call(rbind.data.frame,sapply(ind_results, function(x) { x["SN"] }))))

p_list = list(SNPs="number.of.SNPs.",
              Indels="number.of.indels.",
              Number_Of_MultiAllelic_sites="number.of.multiallelic.sites.",
              Number_Of_MultiAllelic_SNP_Sites="number.of.multiallelic.SNP.sites.")

plots <- lapply(names(p_list), function(x) { 
  
    ggplot(snp_indel[snp_indel$variable==p_list[x],], 
         aes(x=filename, y=value, fill=filename)) +
    geom_bar(stat="identity") +
    geom_text( aes(label=value,y=value, vjust=-1)) +
    facet_grid(. ~ variable) +
    theme(legend.position="bottom") +
    labs(y=x, x="Call Set",title=x) +
    scale_y_continuous(breaks=number_ticks(10), labels=scientific_format())
  
})

#-------------------------#
# Ind. Sample Concordance #
#-------------------------#

con_comb <- do.call("rbind", concordance_results)
con_comb$Concordance <- as.numeric(con_comb$Concordance)

con_comb <- filter(con_comb, Query == Sample) %.%
  group_by(Comparison) %.% 
  mutate(sample_avg_conc=mean(Concordance)) %.%
  arrange(sample_avg_conc)


# Plot concordance of multiple concordance results
ggplot(con_comb) +
    geom_line(position="identity",aes(y=Concordance, x=con_comb$Sample , color=Comparison, group=Comparison)) +
    labs(title="Average Discordance", x="Sample Name", y="Average Discordance (%)") +
    theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(file=sprintf("../reports/quality/%s.png", args[1] ))

#----------------------------#
# Pairwise Concordance Chart #
#----------------------------#

con_chart <- function(record) {
  # Fix data types
  record$Concordance <- as.numeric(record$Concordance)
  
  ggplot() +
    geom_tile(data=record, aes(x=Query, y=Sample, fill=Concordance), colour = "white") + 
    geom_tile(data=record, aes(x=Query, y=Sample, fill=Concordance), colour = "white") + 
    #geom_tile(data=SM_set, aes(x=Sample, y=Sample, fill=Average.Discordance.Number.of.sites), colour = "white") +
    labs(title=sprintf("%s",record[1,'Comparison'])) +
    scale_fill_gradient(low="white", high="red", space="Lab") +
    theme(legend.position="bottom", legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1))
}


# Generate Concordance Charts and place them within a frame.
png("concordance.png", width = 1200, height=round(length(pairs) / 2)*600, units = "px")
do.call(grid.arrange, c(lapply(concordance_results, function(x) { con_chart(x)  }), ncol=2))
dev.off()
