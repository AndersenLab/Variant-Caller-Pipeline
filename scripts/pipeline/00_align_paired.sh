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

#SBATCH --job-name=bwa

#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=16384


#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/fq

BAM_DIR="../bam/"
PICARD="../../tools/"

##################
# Error Checking #
##################


# Check to see if bam exists; if yes, abort
if [ -f "$bam_name" ]
then
	echo "Bam File exists; aborting."
	exit 0
fi

# Check that both files exist.
if [ ! -f "$1" ]
then
	echo "$1 does not exist; aborting."
	exit 0
fi

if [ ! -f "$2" ]
then
	echo "$2 does not exist; aborting."
	exit 0
fi

# Check that a reference is specified.
if [ -z "$3" ]
then
	echo "No reference defined; aborting."
	exit 0
fi

# Double check checksums
f1_md5=`md5sum $1 |  awk '{print substr($0,0,5)}'`
f2_md5=`md5sum $2 |  awk '{print substr($0,0,5)}'`

echo "$1 md5 $f1_md5"
echo "$2 md5 $f2_md5"

if [ "$1" -ne "${f1[3]}" ]; then
	echo "$1; checksum mismatch."
	exit 0
fi

if [ "$2" -ne "${f2[3]}" ]; then
	echo "$2; checksum mismatch."
	exit 0
fi

###############
# Get FASTQ's #
###############

# Convert f1 and f2 to arrays
f1=(${1//-/ })
f2=(${2//-/ })

echo $1
echo $2

# Set reference
ref=${3##*/}

bam_name=${BAM_DIR}${f1[0]}-${f1[1]}-${f1[2]}-${ref%.fa}-${f1[3]}-${f2[3]}.bam

bwa mem -t 8 ../reference/${3} ${1} ${2} | samtools view -b -S -h -  > ${BAM_DIR}${f1[3]}-${f2[3]}.tmp.bam

## Replace Header & Sort Simultaneously!
java	-jar	${PICARD}AddOrReplaceReadGroups.jar	\
I=${BAM_DIR}${f1[3]}-${f2[3]}.tmp.bam	\
O=${BAM_DIR}${f1[3]}-${f2[3]}.tmp.sorted.bam	\
RGID=`basename ${bam_name}`	RGLB=${f1[1]}	\
RGPL=ILLUMINA RGPU=${f1[1]} \
RGSM=${f1[2]}	RGCN=BGI \
SO=coordinate

## Mark Duplicates
java -jar ${PICARD}MarkDuplicates.jar \
I=${BAM_DIR}${f1[3]}-${f2[3]}.tmp.sorted.bam \
O=${bam_name} \
M=${bam_name}.duplicate_report.txt \
VALIDATION_STRINGENCY=SILENT \
REMOVE_DUPLICATES=false

# Remove temp
rm ${BAM_DIR}${f1[3]}-${f2[3]}.tmp.bam ${BAM_DIR}${f1[3]}-${f2[3]}.tmp.sorted.bam


# Index Bams
samtools index $bam_name

