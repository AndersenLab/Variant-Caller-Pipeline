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

debug = T

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


#=======================#
# Get VCF Results       #
#=======================#

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

#=============================#
# Get Concordance Matrix      #
#=============================#

get_concordance_matrix <- function(f1, f2) {
  # Generate a concordance Matrix
  tmp <- tempfile()
  
  # Get sample names from each file
  f1_samples <- unlist(str_split(system(sprintf("bcftools view -h %s| grep '#CHROM' | cut -f '10-'",f1), intern=T),'\t'))
  f2_samples <- unlist(str_split(system(sprintf("bcftools view -h %s| grep '#CHROM' | cut -f '10-'",f2), intern=T),'\t'))
  
  # Get intersection of samples
  joint_samples <- union(f1_samples, f2_samples)
  # Retrieve Concordance Data
  
  # Local Version [MAC]
  if (Sys.info()["sysname"] == "Darwin") {
  r <- read.csv2(pipe(sprintf("echo %s | xargs -t -I \'{}\' -n 1 -P 20 sh -c \"bcftools gtcheck -s {} -G 1 -g %s %s | sed \'s/$/\t{}/\'\" | egrep -v \'#\' ", paste0(joint_samples, collapse=" "), f1, f2)), as.is=T, sep='\t',header=F, col.names= c("CN","Discordance_total", "Discordance_per_site", "Number_of_Sites", "Sample", "Sample_ID", "Query"))
  } else {
  # Cluster Version
  r <- read.csv2(pipe(sprintf("echo -n %s | xargs -d \' \' -t -I \'{}\' -n 1 -P 20 sh -c \"bcftools gtcheck -s {} -G 1 -g %s %s | sed \'s/$/\t{}/\'\" | egrep -v \'#\' ", paste0(joint_samples, collapse=" "), f1, f2)), as.is=T, sep='\t',header=F, col.names= c("CN","Discordance_total", "Discordance_per_site", "Number_of_Sites", "Sample", "Sample_ID", "Query"))
  }
  
  # Fix numbers
  r$Discordance_total <- as.integer(r$Discordance_total)
  r$Discordance_per_site <- NULL
  r$Comparison <- sprintf('%s__%s', repl_multiple(f1, replace_words), repl_multiple(f2, replace_words))

  # Add in Concordance Rate
  r$Concordance <- as.numeric((r$Number_of_Sites - r$Discordance_total) / r$Number_of_Sites)
  r
}

concordance_chart <- function(record, union=F) {
  # Plots concordance. 

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
    #geom_tile(data=SM_set, aes(x=Sample, y=Sample, fill=Average.Discordance.Number.of.sites), colour = "white") +
    labs(title=plot_title, y=str_split(record[1,'Comparison'],"__")[[1]][1], x=str_split(record[1,'Comparison'],"__")[[1]][2]) +
    scale_fill_gradient(low="#FFE900", high="#0092FF", space="Lab", limits=c(0.95,1)) +
    theme(legend.position="bottom", legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, ))
  ggsave(filename = paste0(results_dir, plot_title, ".png", collapse=""), height=20, width=20)
}


args <- sprintf("%s", commandArgs(trailingOnly = TRUE))
setwd("../../data/vcf/")
args <- c("radseq.vcf.gz","test5.20.bcf","test5.30.bcf","test5.40.bcf","test5.50.bcf")

results_dir <- sprintf("../../results/%s/", args[1])
dir.create(results_dir)

pairs <- list()
# args <- c("andersen08_radseq.ws225.vcf.gz", "mmp.vcf.gz")
for (i in 1:(length(args)-1)) {
  pairs<-append(pairs,list(c(args[1], args[i+1])))
}

# Generate Data
vcf_stats <- lapply(args, function(x) { get_vcf_stats(x) })

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

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

lapply( names(snp_indel)[names(snp_indel) != "id"] , function(i) {
  ggplot(snp_indel) + 
    geom_bar(stat="identity", aes_string(x="id", y= i, fill="id")) +
    geom_text( aes_string(x="id",y=i, label= i , vjust=2), color="white", size=12) +
    theme(legend.position="none", legend.title=element_blank()) +
    labs(y=Cap(i), x="Call Set",title=Cap(i) ) +
    scale_y_continuous(breaks=number_ticks(10), labels=comma_format()) + 
    theme(legend.position="none", plot.margin = unit(c(2,2,2,2), "cm")) +
    theme(axis.title.x=element_text(vjust=-2)) +
    theme(axis.title.y=element_text(vjust=4)) +
    theme(plot.title=element_text(size=18, vjust=3)) +
    scale_fill_manual(values=cbPalette) +
    ggsave(filename = paste0(results_dir, i, ".png", collapse=""), plot=last_plot(), height=12)
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

# Individual Sample Concordance
con_comb <- do.call("rbind", concordance_results)

con_comb <- filter(con_comb, Query == Sample) %.%
  group_by(Comparison) %.% 
  mutate(sample_avg_conc=mean(Concordance)) %.%
  ungroup() %.%
  group_by(Sample) %.%
  mutate(sample_count=length(Sample)) %.%
  arrange(Comparison)

# Plot concordance of individual strains
ggplot(con_comb) +
    geom_line(position="identity",aes(y=Concordance, x=con_comb$Sample , color=Comparison, group=Comparison)) + 
    geom_point(position="identity",aes(y=Concordance, x=con_comb$Sample , color=Comparison, group=Comparison)) +
    labs(title="Concordance by Strain", x="Sample Name", y="Concordance (%)") +
    theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1))  +
    facet_grid(. ~ Comparison, drop=T, space="free", scales = "free")

ggsave(filename = paste0(results_dir, "individual_concordance.png"), width=length(last_plot()$data$Query)/7)

## Pairwise Concordance Grid
lapply(concordance_results, function(x) { concordance_chart(x) })
lapply(concordance_results, function(x) { concordance_chart(x, union=T) })


## Plot Concordance by file; compared with the original.
ggplot(concordance_results)

## Save Data Again
save(list = ls(all = TRUE), file= paste0(results_dir, "data.Rdata"))

