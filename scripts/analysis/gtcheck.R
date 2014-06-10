library(ggplot2)
library(gridExtra)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

args <- c("03_RET7a.txt.Q10.vcf.gz","03_RET7a.txt.Q10.vcf.gz","03_RET7a.txt.Q20.vcf.gz","03_RET7a.txt.Q30.vcf.gz","03_RET7a.txt.Q40.vcf.gz")
pairs <- list()
for (i in 1:(length(args)-1)) {
  pairs<-append(pairs,list(c(args[i], args[i+1])))
}

#------------------------------------------#
# Define Functions for generating bcfstats #
#------------------------------------------#

# Function used for plotting later.
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# Pipe in gtcheck data from bcfstats
import_table <- function(t_name, f) { 
  # This function can import data from bcftools stats function
  t <- read.csv2(pipe(sprintf("egrep '^(# \\[1\\])?%s[^,]' %s", t_name, f)), as.is=T, sep='\t',header=T, check.names=FALSE)
  names(t) <- lapply(names(t), function(x) make.names(gsub("\\[[1-9]+\\]","",x)))
  t[,c(-1)]
}


vcf_results <- function(p1,p2) {
  
  
}

get_ind_results <- function(f) {
  # This function retrieves individual data for each VCF for plotting and comparison.
  tmp <- "tmp.txt"
  system(sprintf('bcftools stats -c both -s -  %s > %s', f, tmp))
  tmp_SN  <- import_table("SN", tmp)
  QUAL <- import_table("QUAL", tmp)
  PSC <- import_table("PSC", tmp)
  
  # Fix summary numbers
  SN <- as.list(tmp_SN[[3]])
  names(SN) <- make.names(tmp_SN[[2]])
  print(SN)
  list(SN=SN, QUAL=QUAL, PSC=PSC)
}

# Generate concordance numbers
get_pair_results <- function(f1, f2) {
  # This function retrieves pairwise stats using bcftools for plotting and comparison.
  tmp <- sprintf('%s__%s.txt', sub(".txt","",sub(".vcf.gz","",f1)), sub(".txt","",sub(".vcf.gz","",f2)))
  system(sprintf('bcftools gtcheck  -G -s  %s %s > %s', f1, f2, tmp))
  SM <- import_table("SM", tmp)
  CN <- import_table("CN", tmp)
  SM$comparison <- sprintf('%s__%s', sub(".txt","",sub(".vcf.gz","",f1)), sub(".txt","",sub(".vcf.gz","",f2)))
  
  SM$Average.Discordance.Number.of.sites <- as.numeric(SM$Average.Discordance.Number.of.sites)*100
  
  list(SM=SM,CN=CN)
}

#---------------#
# Generate Data #
#---------------#

ind_results <- lapply(args, function(x) { get_ind_results(x) } )

pair_results <- lapply(pairs, function(x) { get_pair_results(x[1],x[2]) })

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

#-------------------#
# Concordance Chart #
#-------------------#

# Generate Concordance dataframe
SM_frame <- do.call(rbind.data.frame,sapply(pair_results, function(x) { x["SM"] }))

# Generate strain average for sorting.
SM_frame <- group_by(SM_frame, Sample) %.% 
  mutate(sample_avg_depth=mean(Average.Discordance.Number.of.sites)) %.%
  arrange(sample_avg_depth)

# Plot concordance of multiple concordance results
ggplot(SM_frame) +
    geom_line(position="identity",aes(y=Average.Discordance.Number.of.sites, x=reorder(SM_frame$Sample,sample_avg_depth) , color=comparison, group=comparison)) +
    labs(title="Average Discordance", x="Sample Name", y="Average Discordance (%)") +
    theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(breaks = seq(0,100,10))

ggsave(file=sprintf("../reports/quality/%s.png", args[1] ))

#-------------------#
# Concordance Chart #
#-------------------#

con_chart <- function(record) {
  CN_set <- record$CN
  SM_set <- record$SM
  CN_set$Disc <- (CN_set$Discordance/CN_set$Number.of.sites)
  ggplot() +
    geom_tile(data=CN_set, aes(x=Sample.i, y=Sample.j, fill=Disc), colour = "white") + 
    geom_tile(data=CN_set, aes(x=Sample.j, y=Sample.i, fill=Disc), colour = "white") + 
    #geom_tile(data=SM_set, aes(x=Sample, y=Sample, fill=Average.Discordance.Number.of.sites), colour = "white") +
    labs(title=sprintf("%s",SM_set[1,'comparison'])) +
    scale_fill_gradient(low="white", high="red", space="Lab") +
    theme(legend.position="bottom", legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1))
}


# Generate Concordance Charts and place them within a frame.
png("concordance.png", width = 1200, height=round(length(pairs) / 2)*600, units = "px")
do.call(grid.arrange, c(lapply(pair_results, function(x) { con_chart(x)  }), ncol=2))
dev.off()
