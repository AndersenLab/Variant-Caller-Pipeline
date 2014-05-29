# Generate Strain Account Report
library(ggplot2)
library(dplyr)

lib_strain <- read.table("../../data/ancillary/library_strain.txt", header=T)
seq <- read.table("../../data/ancillary/bam_sets/00_all_bams.txt",sep="-", col.names=c("Run","Library","Strain","hash1","hash2"))

# Tabulate

Exp_n <- group_by(lib_strain, Run) %>%
                         mutate(run_count_exp = n()) %>%
                         ungroup() %>%
                         group_by(Run, Library) %>%
                         mutate(lib_count_exp = n()) %>%
                         ungroup() %>%
                         group_by(Strain) %>%
                         mutate(strain_count_exp = n())

Actual_n <- group_by(seq, Run) %>%
                        mutate(run_count_act = n()) %>%
                        ungroup() %>%
                        group_by(Run, Library) %>%
                        mutate(lib_count_act = n()) %>%
                        ungroup() %>%
                        group_by(Strain) %>%
                        mutate(strain_count_act = n())

sum <- unique(left_join(Exp_n, Actual_n, by=c("Run","Library","Strain")))
sum <- select(sum, -contains("hash") )

# Display missing strains

sum[is.na(sum$run_count_act),]

#===========#
# Summarize #
#===========#
# Strains/Lib
r <- unique(sum[!is.na(sum$run_count_act),c("Run", "run_count_exp", "run_count_act")])
excel(r)
# Lib
r <- unique(sum[!is.na(sum$run_count_act),c("Run", "Library", "lib_count_exp", "lib_count_act")])
excel(r)
