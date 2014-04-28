#!bin/bash

#########
# Usage #
#########
#
# Pass 3 parameters - 
# (1) - Search Pattern
# (1) - "subset-on or subset-off"
# (2) - A directory/pattern to search for bam files. For example, if you only wanted to capture BGI2 -
# fasta/*BGI2*.bam would work.


# Combine bam files, sort, re-index, and subset.

for file in `ls $1`
do
	echo $file
	samtools view -b $file chrI:1-100000 > chr1_1st_MB_$file
done