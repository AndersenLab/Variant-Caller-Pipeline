library(VariantAnnotation)

args <- commandArgs(trailingOnly = TRUE)

args <- "test.txt.Q20.vcf.gz"

print(args)

import_table <- function(t_name) { 
  # This function can import data from bcftools stats function
  t <- read.csv2(pipe(sprintf("egrep '^(# )?%s[^,]' %s", t_name, tmp)), sep='\t',header=T, check.names=FALSE)
  names(t) <- lapply(names(t), function(x) gsub("\\[[1-9]+\\]","",x))
  t <- t[,c(-1)]
  t
}

load_stats <- function(vcf) {
# Loads the stats for a given vcf file
tmp <-tempfile()[1]
system(sprintf("bcftools stats %s > %s", vcf, tmp))

# Import tables
AF  <- import_table("AF")      # Stats by Non-Reference Allele
SiS <- import_table("SiS")     # Singleton Stats
ST  <- import_table("ST")      # Substitution Types
tmp_SN  <- import_table("SN")      # Summary Numbers
QUAL<- import_table("QUAL")    # Stats by Quality
IDD <- import_table("IDD")     # InDel Distribution

# Fix summary numbers
SN <- as.list(tmp_SN[[3]])
names(SN) <- tmp_SN[[2]]

# Return list of results.
list(AF=AF,SiS=SiS,ST=ST,SN=SN,QUAL=QUAL,IDD=IDD)
}

q <- load_stats('test.txt.Q40.vcf.gz')


