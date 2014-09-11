#!/usr/bin/python

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


#SBATCH --job-name=mpile

#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mem=32678

#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/bam
from datetime import datetime
startTime = datetime.now()
import os, sys
import os.path

outfile = sys.argv[1]

# Filter depth
fasta_list = open("../ancillary/bam_sets/" + outfile,'r').read().split('\n')

# Initiate Log File
if os.path.isfile("../log2.txt") == False:
	f = open("../log2.txt",'w')
	f.write("jobid\tnode\tsubmission\tdate\tTimeToComplete\tcommand\n")
else:
	f = open('../log2.txt','a+')

f.write("%s\t%s\t%s\t%s\tTimeToComplete_%s\t%s\n" % (os.environ['SLURM_JOB_ID'], os.environ['SLURM_JOB_NODELIST'], outfile, startTime, os.environ['SLURM_JOB_ID'], sys.argv))

f.flush()

# Pull out the chromosomes to so we can split and parallelize
#com = "samtools view -H %s | grep '\@SQ' | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 7 sh -c 'samtools mpileup -r chrIII:1-10000000 -uf ../reference/ce10/ce10.fa -r {} %s | bcftools view -vcg - > ../vcf/tmp.%s.{}.vcf'" % (fasta_list[0], ' '.join(fasta_list),  outfile)
os.chdir("../bam")

# Split~ Chr 3 + Chr 5 for testing.
os.system("cat ../ancillary/chr_ranges_full.txt | xargs -I {} -n 1 -P 12 --verbose sh -c 'samtools mpileup -t DP,DV,DP4,SP -g -D -f ../reference/ce10/ce10.fa -r {} %s | bcftools call --format-fields GQ,GP -m -v > ../vcf/raw.%s.{}.bcf'" % (' '.join(fasta_list),  outfile))
os.system("cat ../ancillary/chr_ranges_full.txt | bcftools concat -O b `ls -v ../vcf/raw.%s.*.bcf` > ../vcf/%s.bcf" % (outfile , outfile.replace(".txt","")))

# Index
os.system("bcftools index -f ../vcf/%s.bcf" % outfile.replace(".txt",""))

# Remove temporary files
os.system("rm -f ../vcf/raw.%s.*" % outfile)

# Replace file with time taken
os.system("sed --in-place='.bak' -e 's/TimeToComplete_%s/%s/' ../log2.txt" % (os.environ['SLURM_JOB_ID'], str(datetime.now()-startTime)))
