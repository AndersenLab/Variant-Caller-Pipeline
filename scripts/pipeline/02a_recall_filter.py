#!/usr/bin/python

#=============#
# Description #
#=============#
#
# This script will submit a job that:
#
# * Recalls variants from the individual pipeline.

#SBATCH --job-name=recall

#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=16384


#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/bam


#=================#
# Process Fastq's #
#=================#

# First - check whether the files in question exist.
import os, glob, sys, subprocess, re, gzip, math, tempfile

sys.path.append("/lscr2/andersenlab/dec211/scripts/pipeline/")
#sys.path.append("/exports/people/andersenlab/dec211/python_modules/peewee/")
#sys.path.remove("/usr/local/lib/python2.6/dist-packages/peewee-2.3.2-py2.6.egg")

print sys.path

from seq_utility_functions import *
from conn import *
from datetime import datetime
from peewee import *


#=======#
# Setup #
#=======#

# Determine system type:
system_type =  os.uname()[0]

# Setup Scripts Dir
scripts_dir = "../../scripts/pipeline/"

# Get Arguments
strain = sys.argv[1]
print strain

#=========#
# Options #
#=========#

# Determine system type:
system_type =  os.uname()[0]
reference = "ce10"
md5_system = {"Linux" : "md5sum", "Darwin":"md5"} # Needed to make md5 work locally and on cluster.
system_cores = {"Linux": 8, "Darwin": 4} # For use with BWA


#========#
# Recall # Incorporate full list of available variants 
#========#

# Generate and store bcf stats
def gen_bcf_stats(strain, bcf, description, Sub_Attribute):
	# Generate stats - Following Het Polarization, Overwrites previous stat file.
	tmp =  tempfile.NamedTemporaryFile().name
	os.system("bcftools stats --samples - %s > %s" % (bcf, tmp))

	# Get variant statistics.
	for i in import_stats("%s" % (tmp)):
		save_eav(strain, i[0], i[1], Entity_Group = "BCF Statistics", Sub_Entity = "%s.%s.bcf" % (strain, description), Sub_Attribute = Sub_Attribute, Tool="bcftools")


# Fetch contigs
contigs = {}
for x in file("../reference/%s/%s.fa.fai" % (reference, reference), 'r').read().split("\n")[:-1]:
	contigs[x.split("\t")[0]] = int(x.split("\t")[1])


command = "parallel --verbose --gnu 'samtools mpileup -t DP,DV,DP4,SP -g -f ../reference/%s/%s.fa -r {} %s.bam | bcftools call -T %s -O b --format-fields GQ,GP -c -v > ../group_bcf/raw.%s.{}.group.bcf' ::: %s" % (reference, reference, strain, '../individual_bcf/complete_variant_set.vcf.gz', strain, ' '.join(contigs.keys()))
os.system(command)
os.system("echo %s | bcftools concat -O b `ls -v ../group_bcf/raw.%s.*.group.bcf` > ../group_bcf/%s.group.bcf" % (contigs.keys(), strain , strain.replace(".txt","")))

# Index
os.system("bcftools index -f ../group_bcf/%s.group.bcf" % strain)

# Remove temporary files
#os.system("rm -f ../group_bcf/raw.%s.*" % strain)


gen_bcf_stats(strain, "../group_bcf/%s.group.bcf" %  strain, "Group", "No Filter")

#========#
# Filter # Variants
#========#
