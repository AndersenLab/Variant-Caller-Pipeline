#!/bin/bash

#=============#
# Description #
#=============#
#
# This script will submit a job that:
#
# * Call variants using mpileup on a specified pattern of bams.
#
# PARAMETERS:
#
# (1) - Pass a pattern to use for bams.


#SBATCH --job-name=depth

#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mem=32768
#SBATCH --mem-per-cpu=2730

#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/bam

echo -e "bam\tavg_depth\tbase_coverage" > ../data/ancillary/bam_depth.txt
ls *.bam | xargs -I {} -n 1 -P 12 sh -c "samtools depth {} | awk -v r={} '{sum+=\$3;cnt++}END{print r\" \"sum/cnt\" \"sum}'  >> ../ancillary/bam_depth.txt" 

ls ../ancillary/bam_sets/ | xargs -I {} -n 1 -P 12 sh -c "samtools depth  -f ../ancillary/bam_sets/{} | awk -v r={} '{sum+=\$3;cnt++}END{print r\" \"sum/cnt\" \"sum}'  >> ../ancillary/bam_depth.txt" 
