#!/usr/bin/python

###############
# Description #
###############
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


outfile = os.path.basename(sys.argv[1])
try:
    options = os.path.basename(sys.argv[2])
except:
    options = []

job_id = os.environ['SLURM_JOB_ID']
os.system("echo \"Job Id: %s\"" % job_id)
os.system("echo \"Set:%s\"" % outfile)

# Depth
# if options.find('d'):
# Calculate the average depth of the bam file.
# avg_depth = os.popen("samtools depth BGI2-RET4-QX1212-11d58-2bc7c.bam  | awk '{sum+=$3;cnt++}END{print sum/cnt}'")


# Filter depth
fasta_list = open("../ancillary/bam_sets/" + outfile,'r').read().split('\n')


print "Job Id: %s" % job_id
print
print "Command: %s" % (sys.argv)
print fasta_list
print

f = file('../ancillary/log_set.txt', 'a+')
f.write("%s - %s\n\n" % (job_id, outfile))

# Pull out the chromosomes to so we can split and parallelize
# com = "samtools view -H %s | grep '\@SQ' | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 7 sh -c 'samtools mpileup -r chrIII:1-10000000 -uf ../reference/ce10/ce10.fa -r {} %s | bcftools view -vcg - > ../vcf/tmp.%s.{}.vcf'" % (fasta_list[0], ' '.join(fasta_list),  outfile)
os.chdir("../bam")

# Split~ Chr 3 + Chr 5 for testing.
os.system("cat ../ancillary/chr_ranges.txt | xargs -I {} -n 1 -P 12 --verbose sh -c 'samtools mpileup -d 100 -D -S -g -f ../reference/ce10/ce10.fa -r {} %s | bcftools view -bvcg - > ../vcf/raw.%s.{}.bcf'" % (' '.join(fasta_list),  outfile))
# awk '{ if ($6>=%s ||  $1 ~ /^#/) print}' - |
os.system("cat ../ancillary/chr_ranges.txt | bcftools cat `ls -v ../vcf/raw.%s.*.bcf` | bcftools view - | vcf-sort > ../vcf/%s.vcf" % (outfile, outfile))
# os.system("bcftools index ../vcf/%s.bcf > ../vcf/%s.bcf" % (outfile, outfile))

# Process VCF
os.system("bgzip -f ../vcf/%s.vcf" % (outfile))
os.system("tabix -f ../vcf/%s.vcf.gz" % (outfile))

# Remove temporary files
os.system("rm -f ../vcf/raw.%s.*" % outfile)

# Print total time
f.write("Time to complete: %s\n\n" % str(datetime.now()-startTime))
