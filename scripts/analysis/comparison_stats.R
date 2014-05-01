library(VariantAnnotation)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

# A function for cleaning variable names; https://github.com/johnmyleswhite/ProjectTemplate/
clean.variable.name <- function(variable.name)
{
  variable.name <- gsub('^[^a-zA-Z0-9]+', '', variable.name, perl = TRUE)
  variable.name <- gsub('[^a-zA-Z0-9]+$', '', variable.name, perl = TRUE)
  variable.name <- gsub('_+', '.', variable.name, perl = TRUE)
  variable.name <- gsub('-+', '.', variable.name, perl = TRUE)
  variable.name <- gsub('\\s+', '.', variable.name, perl = TRUE)
  variable.name <- gsub('\\.+', '.', variable.name, perl = TRUE)
  variable.name <- gsub('[\\\\/]+', '.', variable.name, perl = TRUE)
  variable.name <- make.names(variable.name)
  return(variable.name)
}

#----------------------------#
# Import table from bcfstats #
#----------------------------#

import_table <- function(t_name, tmp) { 
  # This function can import data from bcftools stats function
  t <- read.csv2(pipe(sprintf("egrep '^(# )?%s[^,]' %s", t_name, tmp)), sep='\t',header=T, check.names=FALSE)
  names(t) <- lapply(names(t), function(x) gsub("\\[[1-9]+\\]","",x))
  t <- t[,c(-1)]
  t
}


#--------------------------------------#
# Loads the stats for a given vcf file #
#--------------------------------------#

load_stats <- function(vcf) {
# Loads the stats for a given vcf file
tmp <-tempfile()[1]
system(sprintf("bcftools stats %s > %s", vcf, tmp))

# Import tables
AF  <- import_table("AF", tmp)      # Stats by Non-Reference Allele
SiS <- import_table("SiS", tmp)     # Singleton Stats
ST  <- import_table("ST", tmp)      # Substitution Types
tmp_SN  <- import_table("SN", tmp)  # Summary Numbers
QUAL<- import_table("QUAL", tmp)    # Stats by Quality
IDD <- import_table("IDD", tmp)     # InDel Distribution

# Fix summary numbers
SN <- as.list(tmp_SN[[3]])
names(SN) <- clean.variable.name(tmp_SN[[2]])
# Return list of results.
list(AF=AF,SiS=SiS,ST=ST,SN=SN,QUAL=QUAL,IDD=IDD)
}


# Load multiple vcfs
load_vcfs <- function(search_string) {
  r = list()
  for(var in dir(pattern=sprintf(".*%s.*.vcf.gz$",search_string))) {
    r[[var]] <-load_stats(var)
  }
  r
}

# Function for pulling together variables from 'SN'
pull_var <- function(r, col, var) {
  sapply(r, function(x) {
    x[[col]][[var]]
  })
}

pull_SiS_var <- function(r, var) {
  sapply(r, function(x) {
    x$SiS[[var]]
  })
}


#------#
# Plot #
#------#

# Summary Stats

r <- load_vcfs(args[1])
print(getwd())
# Number of Singletons
x <- as.data.frame(
     cbind(Singletons=pull_var(r,"SiS","number of SNPs")))
x <- cbind(id=rownames(x), x)

c <- ggplot(x, aes(x=id,y=Singletons)) +
    geom_bar(stat="identity") +
    ggtitle(sprintf("%s", args[1] )) +
    stat_bin(geom="text", aes(label=Singletons, vjust=-1))

ggsave(file=sprintf("%s.png", args[2] ))

