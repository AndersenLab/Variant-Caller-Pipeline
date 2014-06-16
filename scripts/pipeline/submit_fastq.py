#!/usr/bin/python

# This script assemble the list of fastqs, identify which have yet to be aligned, and align them (producing BAMs)

# Example Usage:
# sh submit_fastq.sh 

# Prepare necessary libraries

import glob
import os
import sys


#os.system("prepare bwa")
#os.system("prepare samtools")

orig_dir = os.getcwd()

# Change to appropriate directory
os.chdir("../../data/fq")
fasta_list = glob.glob("*.fq.gz")


os.chdir("../bam")

align_fqs = []
fasta_sets = zip(sorted([x for x in fasta_list if "-1.fq.gz" in x]), sorted([x for x in fasta_list if "-2.fq.gz" in x]))
# Check that all fastqs intact, and construct bam list.
for ele in fasta_sets:
	# Check that Run, Library, Sample/Strain are all identicle.
	if (ele[0].split("-")[0:3] == ele[1].split("-")[0:3]) == False:
		sys.exit("Error: Fastqs do not line up.")
	# Check if bam exists
	bam_name = "-".join(ele[0].split("-")[0:4] + [ele[1].split("-")[3]]) + ".bam"
	if os.path.isfile(bam_name) == False:
		align_fqs.append(ele)

# Output to user
if len(align_fqs) > 0:
	print "(%s/%s) Not Aligned" % (len(align_fqs),len(fasta_sets))
	print "Aligning new fastqs:"
	print '\n'.join(["%-50s %-50s" % (x[0], x[1]) for x in align_fqs])
else:
	print "All fastqs are aligned and in Bam format."

# Return to original directory
os.chdir(orig_dir)

# Submit align commands
for f in align_fqs:
	os.system("sbatch 00_align_paired.sh %s %s" % (f[0], f[1]))