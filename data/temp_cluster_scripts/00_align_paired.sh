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
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/fasta

BAM_DIR="../bam/"

###############
# Get FASTQ's #
###############

f1=${1}_1.fq.gz
f2=${1}_2.fq.gz

# Unzip the file; Keep original and force
echo "TIME START aln: $SECONDS"
for fa in ${f1} ${f2}
do
	pwd
	bwa aln -t 8 ../reference/ce10/ce10.fa $fa > ${fa/.fq.gz/.sai}
done

bwa sampe ../reference/ce10/ce10.fa ${f1/.fq.gz/.sai} ${f2/.fq.gz/.sai} ${f1} ${f2} | samtools view -b -S -h - | samtools sort - $BAM_DIR`basename ${f1/_1.fq.gz/}`
# Remove intermediates
rm ${f1/.fq.gz/.sai}
rm ${f2/.fq.gz/.sai}

echo "TIME END aln: $SECONDS"
echo "${BAM_DIR}${f1/.fq.gz/mem}"
echo "TIME START mem: $SECONDS"
bwa mem -t 8 ../reference/ce10/ce10.fa ${f1} ${f2} | samtools view -b -S -h  - | samtools sort - ${BAM_DIR}`basename ${f1/.fq.gz/mem}`
echo "TIME END mem: $SECONDS"

# Index Bams
samtools index $BAM_DIR`basename ${f1/_1.fq.gz/}`.bam
samtools index ${BAM_DIR}`basename ${f1/.fq.gz/mem}`.bam

