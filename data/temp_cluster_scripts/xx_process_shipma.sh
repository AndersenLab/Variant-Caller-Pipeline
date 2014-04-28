#!/bin/bash

#=============#
# Description #
#=============#
#
# This script will submit a job that:
#
# * Aligns, creates bams
# * Sorts
# * Indexes

#SBATCH --job-name=bam

#SBATCH --output=../../log/%j.txt
#SBATCH --error=../../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=16384
#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/bam/shipma_bam

BAM_DIR="../bam/"
PICARD="../../../tools/"

f1=(${1//_/ })
echo ${1}


## Replace Header & Sort Simultaneously!
java	-jar	${PICARD}AddOrReplaceReadGroups.jar	\
I=${1}	\
O=${1/.bam/.tmp.sorted.bam}	\
RGID=${1}	RGLB=${f1[1]}	\
RGPL=ILLUMINA RGPU=${f1[1]} \
RGSM=${f1[2]/.bam/}	RGCN=BGI \
SO=coordinate \
VALIDATION_STRINGENCY=LENIENT

## Mark Duplicates
java -jar ${PICARD}MarkDuplicates.jar \
I=${1/.bam/.tmp.sorted.bam} \
O=${1}.fixed.bam \
M=${1}.duplicate_report.txt \
VALIDATION_STRINGENCY=SILENT \
REMOVE_DUPLICATES=false

# Remove temp
#rm ${1}
rm ${1/.bam/.tmp.sorted.bam}


# Index Bams
samtools index ${1}.fixed.bam

