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


number_ticks <- function(n) {function(limits) pretty(limits, n)}


#=======================#
# Get Pair Results      #
#=======================#
get_pair_results <- function(f1, f2) {
  # This function retrieves individual data for each VCF for plotting and comparison.
  tmp <- tempfile()
  sprintf('bcftools stats -c both -s -  %s %s > %s', f1, f2, tmp)
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

replace_words = c(".vcf.gz",".txt", "../../data/vcf/")

repl_multiple <- function(string, replacements) {
  # Convenience Function for replacing multiple string with nothing.
  for (rep in replacements) {
    string <- sub(rep, "", string)
  }
  string
}

get_con_matrix <- function(f1, f2) {
  # Generate a concordance Matrix
  tmp <- tempfile()
  
  # Get sample names from each file
  f1_samples <- unlist(str_split(system(sprintf("bcftools view -h %s| grep '#CHROM' | cut -f '10-'",f1), intern=T),'\t'))
  f2_samples <- unlist(str_split(system(sprintf("bcftools view -h %s| grep '#CHROM' | cut -f '10-'",f2), intern=T),'\t'))
  
  # Get intersection of samples
  joint_samples <- union(f1_samples, f2_samples)
  
  # Retrieve Concordance Data
  r <- read.csv2(pipe(sprintf("echo %s | xargs -I {} -n 1 -P 50 sh -c \"bcftools gtcheck  -s {} -G 1 -g %s %s | sed \'s/$/\t{}/\'\" | egrep -v '#' ", paste0(joint_samples, collapse=" "), f1, f2)), as.is=T, sep='\t',header=F, col.names= c("CN","Discordance_total", "Discordance_per_site", "Number_of_Sites", "Sample", "Sample_ID", "Query"))
  
  # Fix numbers
  r$Discordance_total <- as.integer(r$Discordance_total)
  r$Discordance_per_site <- NULL
  r$Comparison <- sprintf('%s__%s', repl_multiple(f1, replace_words), repl_multiple(f2, replace_words))
  # Add in Concordance Rate
  r$Concordance <- ((r$Number_of_Sites - r$Discordance_total) / r$Number_of_Sites)
  r
}

con_chart <- function(record) {
  # Fix data types
  record$Concordance <- as.numeric(record$Concordance)
  
  # Use to subset records in common
  # record <- filter(record, Sample %in% Query & Query%in%Sample)
  
  ggplot() +
    geom_tile(data=record, aes(x=Query, y=Sample, fill=Concordance), colour = "white") + 
    geom_tile(data=record, aes(x=Query, y=Sample, fill=Concordance), colour = "white") + 
    #geom_tile(data=SM_set, aes(x=Sample, y=Sample, fill=Average.Discordance.Number.of.sites), colour = "white") +
    labs(title=sprintf("%s",record[1,'Comparison']), y=str_split(record[1,'Comparison'],"__")[[1]][1], x=str_split(record[1,'Comparison'],"__")[[1]][2]) +
    scale_fill_gradient(low="white", high="red", space="Lab") +
    theme(legend.position="bottom", legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, ))
  
}

con_chart_union <- function(record) {
  # Fix data types
  record$Concordance <- as.numeric(record$Concordance)
  
  # Use to subset records in common
  record <- filter(record, Sample %in% Query & Query%in%Sample)
  
  ggplot() +
    geom_tile(data=record, aes(x=Query, y=Sample, fill=Concordance), colour = "white") + 
    geom_tile(data=record, aes(x=Query, y=Sample, fill=Concordance), colour = "white") + 
    #geom_tile(data=SM_set, aes(x=Sample, y=Sample, fill=Average.Discordance.Number.of.sites), colour = "white") +
    labs(title=sprintf("%s",record[1,'Comparison']), y=str_split(record[1,'Comparison'],"__")[[1]][1], x=str_split(record[1,'Comparison'],"__")[[1]][2]) +
    scale_fill_gradient(low="white", high="red", space="Lab") +
    theme(legend.position="bottom", legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1, ))
  
}

