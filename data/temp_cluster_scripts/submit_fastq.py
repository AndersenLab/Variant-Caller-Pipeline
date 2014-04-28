#!/usr/bin/python

# This script will submit a series of fasta files (located within a subfolder) for alignment. 

# Example Usage:
# bash submit_fastq.sh IndexCHKPEI13010003 - All of the files in this directory will be aligned.

# Prepare necessary libraries

import os

# prepare bwa,etc.
os.system("prepare bwa")
os.system("prepare samtools")

for line in  open('../ancillary/fasta_list.txt','r'):
	f = line.replace('\n','').split(',')
	# f = [x for x in f if x.find("QG557")!=-1] # Subset
	if len(f) > 0:
		os.system("sbatch 00_align_paired.sh %s %s" % (f[0], f[1]))