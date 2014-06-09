#!/usr/bin/python

# This script will submit a series of fasta files (located within a subfolder) for alignment. 

# Example Usage:
# bash submit_fastq.sh IndexCHKPEI13010003 - All of the files in this directory will be aligned.

# Prepare necessary libraries

import os
import sys
import subprocess

# prepare bwa,etc.
os.system("prepare bwa")
os.system("prepare samtools")


if len(sys.argv) == 1:
	grep = "."
else:
	grep = sys.argv[1]

file_list = [x.split(',') for x in filter(len,os.popen('grep %s ../../data/ancillary/fastq_list.txt' % grep).read().split("\n")) if x.find("JU1652") != -1]

for f in  file_list:
	os.system("echo %s; echo %s" % (f[0],[f[1]]))
	os.system("sbatch 00_align_paired.sh %s %s" % (f[0], f[1]))