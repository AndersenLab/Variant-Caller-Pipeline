#!/usr/bin/python

#=============#
# Description #
#=============#
#
# This script will submit a job that:
#
# * Merges variants from the individual_bcf folder for use in variant calling.

#SBATCH --job-name=merind

#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=16384


#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/individual_bcf

import os, glob, sys
sys.path.append("/lscr2/andersenlab/dec211/scripts/pipeline/")
from conn import *

#===================#
# Conduct bcf merge #
#===================#
files = glob.glob("*.bcf")

# record member files and count.
save_eav("Individual BCF Merging", "Member bcfs", files, Entity_Group = "BCF Statistics", Tool="bcftools")
save_eav("Individual BCF Merging", "Member bcfs count", len(files), Entity_Group = "BCF Statistics", Tool="bcftools")

os.system("bcftools merge -O z `ls *.single.bcf` > complete_variant_set.vcf.gz")
os.system("bcftools index -f complete_variant_set.bcf")
