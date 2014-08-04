#!/usr/bin/Rscript

#=============#
# Description #
#=============#
#
# This script will submit a job that:
#
# * Compares two vcf files. Requires bcftools 0.2.0-rc8
#
#
#SBATCH --job-name=vcfcomp

#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mem=16384

#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/vcf

library(ggplot2)
library(scales)
library(stringr)
library(dplyr)
library(grid)


# Functions for working with bcftools + VCF
#=======================#
# Import BCFTools Table #
#=======================#
import_table <- function(t_name, f) { 
  # This function can import data from bcftools stats function
  t <- read.csv2(pipe(sprintf("egrep '^(# )%s[^,]|^%s' %s", t_name, t_name, f)), as.is=T, sep='\t',header=T, check.names=FALSE)
  names(t) <- lapply(names(t), function(x) make.names(gsub("\\[[1-9]+\\]","",x)))
  # Remove first column and return table.
  t[,c(-1)]
}


#==========================#
# Ancillary Functions/Defs #
#==========================#

number_ticks <- function(n) {function(limits) pretty(limits, n)}

repl_multiple <- function(string, replacements) {
  # Convenience Function for 
  # replacing multiple string with nothing.
  for (rep in replacements) {
    string <- sub(rep, "", string)
  }
  string
}

replace_words = c(".vcf.gz",".txt", "../../data/vcf/",".bcf")


#=================#
# Fetch Data      #
#=================#

# Get individual VCF Stats
get_vcf_stats <- function(f) {
  # This function retrieves comparison data for
  # each VCF.
  tmp <- tempfile()
  f_name <- repl_multiple(f,replace_words)
  system(sprintf('bcftools stats -s - %s > %s', f, tmp))
  SN  <- import_table("SN", tmp)
  QUAL <- import_table("QUAL", tmp)
  PSC <- import_table("PSC", tmp)
  
  comp <- sprintf('%s', f_name)
  
  # Replace ID Field as is appropriate
  QUAL$id <- f_name
  PSC$id <- f_name
  
  # Fix SN Group, ID
  SN <- reshape(SN, direction="wide", timevar=c("key"), ids=c("id"))
  names(SN) <- make.names(gsub("value.","",gsub(":","", names(SN))))
  SN$id <- f_name
  
  list(SN=SN, QUAL=QUAL, PSC=PSC)
}

# Get concordance results
get_concordance_matrix <- function(f1, f2) {
  # Get sample names from each file
  f1_samples <- unlist(str_split(system(sprintf("bcftools view -h %s| grep '#CHROM' | cut -f '10-'",f1), intern=T),'\t'))
  f2_samples <- unlist(str_split(system(sprintf("bcftools view -h %s| grep '#CHROM' | cut -f '10-'",f2), intern=T),'\t'))
  
  
  # Retrieve Intersection & Absolute Counts.
  tmp <- tempfile()
  system(sprintf('bcftools stats -c both -s -  %s %s > %s', f1, f2, tmp))
  SN  <- import_table("SN", tmp)
  
  SN <- reshape(SN, direction="wide", timevar=c("key"), ids=c("id"))
  names(SN) <- make.names(gsub("value.","",gsub(":","", names(SN))))
  
  
  # Local Version [MAC]
  if (Sys.info()["sysname"] == "Darwin") {
    r <- read.csv2(pipe(sprintf("echo %s | xargs -t -I \'{}\' -n 1 -P 20 sh -c \"bcftools gtcheck -s {} -G 1 -g %s %s | sed \'s/$/\t{}/\'\" | egrep -v \'#\' ", paste0(joint_samples, collapse=" "), f1, f2)), as.is=T, sep='\t',header=F, col.names= c("CN","Discordance_total", "Discordance_per_site", "Number_of_Sites", "Sample", "Sample_ID", "Query"))
  } else {
    save(list = ls(all = TRUE), file= paste0(results_dir, "data2.Rdata"))
    system(sprintf("echo -n %s | xargs -d \' \' -t -I \'{}\' -n 1 -P 20 sh -c \"bcftools gtcheck -s {} -G 1 -g %s %s | sed \'s/$/\t{}/\' | egrep -v \'#\' > '%s{}.CN.txt'\"", paste0(f2_samples, collapse=" "), f1, f2, results_dir))
    r <- lapply(f2_samples, function(x) {
      read.table(sprintf("%s%s.CN.txt", results_dir,x),as.is=T, sep='\t',header=F, col.names= c("CN","Discordance_total", "Discordance_per_site", "Number_of_Sites", "Sample", "Sample_ID", "Query"))
    })
    r <- do.call("rbind", r)
  }
  
  save(list = ls(all = TRUE), file= paste0(results_dir, "data2.Rdata"))
  
  # Fix numbers
  r$Discordance_total <- as.integer(r$Discordance_total)
  r$Number_of_Sites <- as.integer(r$Number_of_Sites)
  r$isec_sites <- SN[,"number.of.SNPs"][3]
  r$f1_sites <- SN[,"number.of.SNPs"][1]
  r$f2_sites <- SN[,"number.of.SNPs"][2]
  r$Discordance_per_site <- NULL
  r$Comparison <- sprintf('%s__%s', repl_multiple(f1, replace_words), repl_multiple(f2, replace_words))
  
  # Add in Concordance Rate
  r$Concordance <- as.numeric((r$Number_of_Sites - r$Discordance_total) / r$Number_of_Sites)
  r$Union_Concordance <- as.numeric(r$Number_of_Sites - r$Discordance_total) / sum(SN[,"number.of.SNPs"])
  SNP_counts <- SN[,"number.of.SNPs"]
  r
}

#========#
# Charts #
#========#

concordance_chart <- function(record, union=F) {
  # Use to subset records in common
  if (union==T) {
    record <- filter(record, Sample %in% Query & Query%in%Sample)
    plot_title <- sprintf("%s_Union",record[1,'Comparison'])
  } else {
    plot_title <- sprintf("%s",record[1,'Comparison'])
  }
  
  ggplot() +
    geom_tile(data=record, aes(x=Query, y=Sample, fill=Concordance), colour = "white") + 
    geom_tile(data=record, aes(x=Query, y=Sample, fill=Concordance), colour = "white") + 
    geom_text(data=record, aes(x=Query, y=Sample, label=format(Concordance*100,digits=3)), colour = "white") + 
    #geom_tile(data=SM_set, aes(x=Sample, y=Sample, fill=Average.Discordance.Number.of.sites), colour = "white") +
    labs(title=plot_title, y=str_split(record[1,'Comparison'],"__")[[1]][1], x=str_split(record[1,'Comparison'],"__")[[1]][2]) +
    scale_fill_gradient(low="#FFE900", high="#0092FF", space="Lab") +
    theme(legend.position="bottom", legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, ))
  ggsave(filename = paste0(results_dir, plot_title, ".png", collapse=""), height=20, width=15+length(args))
}


venn <- function(t="SNP") {
  grid.newpage()
  venn.plot <- draw.pairwise.venn(area1        = sum(complete.cases(jvcf[TYPE == t ,c("REF.x"), with = F])),
                                  area2        = sum(complete.cases(jvcf[TYPE == t,c("REF.y"), with = F])),
                                  cross.area   = sum(complete.cases(jvcf[TYPE == t,c("REF.x","REF.y"), with = F])),
                                  scaled       = T,
                                  category     = str_replace_all(c(f1, f2), c("(.vcf|.gz|.txt|.bcf)"),""),
                                  fill         = c("blue", "red"),
                                  alpha        = 0.3,
                                  lty          = "blank",
                                  cex          = 2,
                                  cat.cex      = 1.5,
                                  cat.pos      = c(0, 0),
                                  cat.dist     = 0.05,
                                  # cat.just     = list(c(-1, -1), c(1, 1)),
                                  ext.pos      = 20,
                                  ext.dist     = -0.05,
                                  ext.length   = 0.85,
                                  ext.line.lwd = 5,
                                  ext.line.lty = "dashed")
  grid.draw(venn.plot)
}


#========#
# Start  #
#========#

args <- sprintf("%s", commandArgs(trailingOnly = TRUE))
setwd("../../data/vcf/")

results_dir <- sprintf("../../results/%s/", args[1])
dir.create(results_dir)

pairs <- list()
# args <- c("andersen08_radseq.ws225.vcf.gz", "mmp.vcf.gz")
for (i in 2:(length(args)-1)) {
  pairs<-append(pairs,list(c(args[2], args[i+1])))
}

print(pairs)

#Generate Data
vcf_stats <- lapply(args[2:length(args)], function(x) { get_vcf_stats(x) })

#-----------------------------------------------#
# Plot Individual vcf results for SNPs + Indels #
#-----------------------------------------------#

# Bring together total SNP counts
snp_indel <- do.call(rbind.data.frame,sapply(vcf_stats, function(x) { x["SN"] }))


# Variables to plot
p_list <- names(snp_indel)[names(snp_indel) != "id"]

# Nice Labels
Cap <- function(x) {
  s <- strsplit(x, " |\\.")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

cbPalette <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#A1887F", "#FFB300", "#0080ff","#0080ff", "#408000"), 5)

lapply( names(snp_indel)[names(snp_indel) != "id"] , function(i) {
  ggplot(snp_indel) + 
    geom_bar(stat="identity", aes_string(x="id", y= i, fill="id")) +
    geom_text( aes_string(x="id",y=i, label= i , vjust=2), color="white", size=9) +
    theme(legend.position="none", legend.title=element_blank()) +
    labs(y=Cap(i), x="Call Set",title=Cap(i) ) +
    scale_y_continuous(breaks=number_ticks(10), labels=comma_format()) + 
    theme(legend.position="none", plot.margin = unit(c(2,2,2,2), "cm")) +
    theme(axis.title.x=element_text(vjust=-2)) +
    theme(axis.title.y=element_text(vjust=4)) +
    theme(plot.title=element_text(size=18, vjust=3)) +
    scale_fill_manual(values=cbPalette) +
    ggsave(filename = paste0(results_dir, i, ".png", collapse=""), plot=last_plot(), height=12, width=15+length(args))
})

# Number of SNPs/Strain
# Plot SNP data by Strain #
PSC_data <- do.call(rbind.data.frame,sapply(vcf_stats, function(x) { x["PSC"] }))
PSC_data <- mutate(PSC_data, total_snps = nRefHom + nNonRefHom + nHets) %.%
  mutate(nonRef_snps = nHets + nNonRefHom)

# Total SNPs
ggplot(PSC_data) + 
  geom_point(aes(x=sample, y= nonRef_snps, group=id, color=id), size=5) +
  labs(title="SNP count/file", x="Strain Name", y="non-reference Allele") +
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values=cbPalette) +
  ggsave(filename = paste0(results_dir, "Total_SNPs.png", collapse=""), width=14)

# Singletons
ggplot(PSC_data) + 
  geom_point(aes(x=sample, y=nSingletons, group=id, color=id), size=5) +
  labs(title="Singletons", x="Strain Name", y="[Singletons] non-reference Allele") +
  theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1))  +
  scale_color_manual(values=cbPalette) +
  ggsave(filename = paste0(results_dir, "Singletons_SNPs.png", collapse=""),  width=14)

save(list = ls(all = TRUE), file= paste0(results_dir, "data.Rdata"))


#---------------------#
# Examine Concordance #
#---------------------#
concordance_results <- lapply(pairs, function(x) { get_concordance_matrix(x[1],x[2]) })

save(list = ls(all = TRUE), file= paste0(results_dir, "data2.Rdata"))

# Individual Sample Concordance
con_comb <- do.call("rbind", concordance_results)

con_comb <- filter(con_comb, Query == Sample) %.%
  group_by(Comparison) %.% 
  mutate(sample_avg_conc=mean(Concordance)) %.%
  ungroup() %.%
  group_by(Sample) %.%
  mutate(sample_count=length(Sample)) %.%
  arrange(Comparison)

con_comb$QUAL <- as.integer(str_match(con_comb$Comparison,"Q([0-9]+)")[,2])
con_comb$Depth <- as.integer(str_match(con_comb$Comparison,"d([0-9]+)")[,2])

# Plot concordance of individual strains
ggplot(con_comb) +
  geom_line(position="identity",aes(y=Concordance, x=con_comb$Sample , color=Comparison, group=Comparison)) + 
  geom_point(position="identity",aes(y=Concordance, x=con_comb$Sample , color=Comparison, group=Comparison)) +
  labs(title="Concordance by Strain", x="Sample Name", y="Concordance (%)") +
  theme(legend.position="right", axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = paste0(results_dir, "individual_concordance.png"), width=10)

ggplot(con_comb) +
  geom_line(position="identity",aes(y=Union_Concordance, x=con_comb$Sample , color=Comparison, group=Comparison)) + 
  geom_point(position="identity",aes(y=Union_Concordance, x=con_comb$Sample , color=Comparison, group=Comparison)) +
  labs(title="Concordance by Strain", x="Sample Name", y="Concordance (%)") +
  theme(legend.position="right", axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = paste0(results_dir, "individual_concordance_union.png"), width=10)



# Fix up df of con_comb.
con_comb <- group_by(con_comb, Comparison, Sample)


fix_labels <- function(x) {
  sub(".bcf","",str_split_fixed(x,"__",2)[,2])
}

  
  # Stratified Concordance
  ggplot(con_comb) +
    geom_line(position="identity",aes(x=Comparison, y=con_comb$Union_Concordance , color=Sample, group=Sample), alpha=0.5) +
    stat_summary(fun.y=mean, mapping = aes(x=con_comb$Comparison, group="Comparison", y = con_comb$Union_Concordance), geom="line", size = 2) +
    labs(title=sprintf("%s vs:Concordance by Strain",args[2]), x="Sample Name", y="Concordance (%)") +
    theme(legend.position="right", legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, )) +
    scale_x_discrete(labels=fix_labels) +
    ggsave(filename = paste0(results_dir, "stratified_concordance_union.png"), width=14)
  
  
  # Plot QUAL
  if (length(is.na(con_comb$QUAL)[is.na(con_comb$QUAL)==F]) > 0) {
    # Plot QUAL if possible.
    qual_set <- con_comb
    qual_set$QUAL[is.na(con_comb$QUAL) == T] <- 0
    ## Plot Qualities within range compared with the the first vcf input.
    ggplot(qual_set) +
      geom_line(position="identity",aes(x=QUAL, y=con_comb$Union_Concordance , color=Sample, group=Sample), alpha=0.5) +
      geom_point(position="identity",aes(x=QUAL, y=con_comb$Union_Concordance , color=Sample, group=Sample), alpha=0.5) +
      stat_summary(fun.y=mean, mapping = aes(x=QUAL, group="Comparison", y = con_comb$Union_Concordance), geom="line", size = 2) +
      labs(title=sprintf("%s vs %s filtered by quality",args[2], args[3]), x="Quality", y="Concordance (%)") +
      theme(legend.position="right", legend.position="top", axis.text.x = element_text(hjust = 1)) +
      scale_x_continuous(breaks=unique(qual_set$QUAL)) +
      ggsave(filename = paste0(results_dir, "stratified_concordance_QUAL_union.png"), width=14)
  }
  
  # Plot DEPTH
  if (length(is.na(con_comb$Depth)[is.na(con_comb$Depth)==F]) > 0) {
    # Plot QUAL if possible.
    depth_set <- con_comb
    depth_set$Depth[is.na(con_comb$Depth) == T] <- NA
    ## Plot Qualities within range compared with the the first vcf input.
    ggplot(depth_set) +
      geom_line(position="identity",aes(x=Depth, y=con_comb$Union_Concordance , color=Sample, group=Sample), alpha=0.5) +
      geom_point(position="identity",aes(x=Depth, y=con_comb$Union_Concordance , color=Sample, group=Sample), alpha=0.5) +
      stat_summary(fun.y=mean, mapping = aes(x=Depth, group="Comparison", y = con_comb$Union_Concordance), geom="line", size = 2) +
      labs(title=sprintf("%s vs %s filtered by quality",args[2], args[3]), x="Depth", y="Concordance (%)") +
      theme(legend.position="right", legend.position="top", axis.text.x = element_text(hjust = 1)) +
      scale_x_continuous(breaks=unique(depth_set$Depth)) +
      ggsave(filename = paste0(results_dir, "stratified_concordance_Depth_union.png"), width=14)
  }
  


## Pairwise Concordance Grid
#lapply(concordance_results, function(x) { concordance_chart(x) })
# Pairwise Concordance Charts
lapply(concordance_results, function(x) { concordance_chart(x, union=T) })

## Save Data Again
save(list = ls(all = TRUE), file= paste0(results_dir, "data.Rdata"))

