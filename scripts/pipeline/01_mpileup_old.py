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
#SBATCH --mem=32768
#SBATCH --mem-per-cpu=2730

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


com = "cat ../ancillary/chr_ranges.txt | xargs -I {} -n 1 -P 12 sh -c 'samtools mpileup -d 10000 -D -S -gu -f ../reference/ce10/ce10.fa -r {} %s | bcftools view -bvcg - > ../vcf/raw.%s.{}.bcf'" % (' '.join(fasta_list),  outfile)

os.system(com)
# awk '{ if ($6>=%s ||  $1 ~ /^#/) print}' - |
os.system("cat ../ancillary/chr_ranges.txt  | bcftools cat `ls -v ../vcf/raw.%s.*.bcf` | bcftools view - | vcf-sort > ../vcf/%s.vcf" % (outfile , outfile))
#os.system("bcftools index ../vcf/%s.bcf > ../vcf/%s.bcf" % (outfile, outfile))

# Process original file
os.system("bgzip -f ../vcf/%s.vcf" % (outfile))
os.system("tabix -f ../vcf/%s.vcf.gz" % (outfile))


# Remove temporary files
os.system("rm -f ../vcf/raw.%s.*" % outfile)

# Print total time
# Replace file with time taken
os.system("sed --in-place='.bak' -e 's/TimeToComplete_%s/%s/' ../log2.txt" % (os.environ['SLURM_JOB_ID'], str(datetime.now()-startTime)))
