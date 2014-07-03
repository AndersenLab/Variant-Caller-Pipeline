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

replace_words = c(".vcf.gz",".txt", "../../data/vcf/")


#=======================#
# Get Pair Results      #
#=======================#

get_pair_results <- function(f1, f2) {
  # This function retrieves comparison data for
  # each VCF.
  tmp <- tempfile()
  system(sprintf('bcftools stats -c both -s -  %s %s > %s', f1, f2, tmp))
  SN  <- import_table("SN", tmp)
  QUAL <- import_table("QUAL", tmp)
  PSC <- import_table("PSC", tmp)
    
  comp <- sprintf('%s__%s', repl_multiple(f1,replace_words), repl_multiple(f2,replace_words))
  id_repl <- list('0'=f1, '1'=f2, '2'=comp)
  
  # Replace ID Field as is appropriate
  SN$key <- tolower(SN$key)
  QUAL$id <- as.character(id_repl[as.character(QUAL$id)])
  PSC$id <- as.character(id_repl[as.character(PSC$id)])
  
  # Fix SN Group
  SN <- reshape(SN, direction="wide", timevar=c("key"), ids=c("id"))
  names(SN) <- make.names(gsub("value.","",gsub(":","", names(SN))))
  
  # Cleanup Names
  SN$name <- repl_multiple(as.character(id_repl[as.character(SN$id)]), replace_words)
  
  # Add union to subtraction (e.g. 2+1->1; 0+2->0); SNP counts don't represent subtraction - instead
  # show them as actual counts.
  SN_vars = c("number.of.snps", "number.of.mnps", "number.of.indels", "number.of.others", "number.of.multiallelic.sites","number.of.multiallelic.snp.sites")
  for (i in SN_vars) {
    SN[SN$id==0,i] <- sum(SN[SN$id==0 | SN$id == 2, i]) # First VCF
    SN[SN$id==1,i] <- sum(SN[SN$id==1 | SN$id == 2, i]) # Second VCF
  }
  
  # Add Comparison Column
  SN$Comparison <- comp
  QUAL$Comparison <- comp
  PSC$Comparison <- comp
  
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
  r <- read.csv2(pipe(sprintf("echo -n %s | xargs -d \' \' -t -I \'{}\' -n 1 -P 20 sh -c \"bcftools gtcheck -s {} -G 1 -g %s %s | sed \'s/$/\t{}/\'\" | egrep -v \'#\' ", paste0(joint_samples, collapse=" "), f1, f2)), as.is=T, sep='\t',header=F, col.names= c("CN","Discordance_total", "Discordance_per_site", "Number_of_Sites", "Sample", "Sample_ID", "Query"))
  
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

results_dir <- sprintf("../../results/%s/", args[1])
dir.create(results_dir)

pairs <- list()
# args <- c("andersen08_radseq.ws225.vcf.gz", "mmp.vcf.gz")
for (i in 2:(length(args)-1)) {
  pairs<-append(pairs,list(c(args[i], args[i+1])))
}

# Generate Data
pair_results <- lapply(pairs, function(x) { get_pair_results(x[1],x[2]) })
concordance_results <- lapply(pairs, function(x) { get_concordance_matrix(x[1],x[2]) })

#-----------------------------------------------#
# Plot Individual vcf results for SNPs + Indels #
#-----------------------------------------------#

# Bring together total SNP counts
snp_indel <- do.call(rbind.data.frame,sapply(pair_results, function(x) { x["SN"] }))

p_list = list(SNPs="number.of.snps",
              Indels="number.of.indels",
              Number_Of_MultiAllelic_sites="number.of.multiallelic.sites",
              Number_Of_MultiAllelic_SNP_Sites="number.of.multiallelic.snp.sites")

lapply(names(p_list), function(i) {
  ggplot(snp_indel) + 
    geom_bar(stat="identity", aes_string(x="name", y=p_list[[i]])) +
    geom_text( aes_string(x="name",y=p_list[[i]], label=p_list[[i]], vjust=2), color="white") +
    theme(legend.position="bottom", legend.title=element_blank()) +
    labs(y=i, x="Call Set",title=i) +
    scale_y_continuous(breaks=number_ticks(10), labels=scientific_format()) + 
    facet_grid(. ~ Comparison, space = "free", scales="free", drop=T ) +
    theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(filename = paste0(results_dir,p_list[[i]], ".png", collapse=""), plot=last_plot(), height=12)
})

# Number of SNPs/Strain
# Plot SNP data by Strain #
PSC_data <- do.call(rbind.data.frame,sapply(pair_results, function(x) { x["PSC"] }))
PSC_data <- mutate(PSC_data, total_snps = nRefHom + nNonRefHom + nHets) %.%
            mutate(nonRef_snps = nHets + nNonRefHom) %.%
            filter(nRefHom > 0) %.%
            filter(id == Comparison) %.%
            arrange(Comparison)

# Total SNPs
ggplot(PSC_data) + 
    geom_point(aes(x=sample, y= nonRef_snps, group=id, color=id), size=5) +
    labs(title="SNP count/file", x="Strain Name", y="non-reference Allele") +
    theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_grid(. ~ Comparison, drop=T, space="free", scales = "free")

    ggsave(filename = paste0(results_dir, "Total_SNPs.png", collapse=""), width=length(last_plot()$data$sample)/7)

# Singletons
ggplot(PSC_data) + 
    geom_point(aes(x=sample, y=nSingletons, group=id, color=id), size=5) +
    labs(title="Singletons", x="Strain Name", y="non-reference Allele") +
    theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1))  +
    facet_grid(. ~ Comparison, drop=T, space="free", scales = "free")

    ggsave(filename = paste0(results_dir, "Singletons_SNPs.png", collapse=""),  width=length(last_plot()$data$sample)/7)

save(list = ls(all = TRUE), file= paste0(results_dir, "data.Rdata"))


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

## Pairwise Concordance
lapply(concordance_results, function(x) { concordance_chart(x) })
lapply(concordance_results, function(x) { concordance_chart(x, union=T) })

## Save Data!
save(list = ls(all = TRUE), file= paste0(results_dir, "data.Rdata"))

