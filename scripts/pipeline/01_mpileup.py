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

os.system("cat ../ancillary/chr3ranges.txt")
#job_id = os.environ['SLURM_JOB_ID']
job_id = "1"


outfile = os.path.basename(sys.argv[1]).replace('[','\[').replace(']','\]')


fasta_list = open("../ancillary/fasta_sets/" + outfile,'r').read().split('\n')


print "Job Id: %s" % job_id
print
print "Command: %s" % (sys.argv)
print fasta_list
print

f = file('../ancillary/log_set.txt','a+')
f.write("%s - %s\n\n" % (job_id, outfile))

# Pull out the chromosomes to so we can split and parallelize
#com = "samtools view -H %s | grep '\@SQ' | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 7 sh -c 'samtools mpileup -r chrIII:1-10000000 -uf ../reference/ce10/ce10.fa -r {} %s | bcftools view -vcg - > ../vcf/tmp.%s.{}.vcf'" % (fasta_list[0], ' '.join(fasta_list),  outfile)
os.chdir("../bam")
# Split~ Chr 3 only, for testing.

com = "cat ../ancillary/chr3ranges.txt | xargs -I {} -n 1 -P 12 sh -c 'samtools mpileup -d 10000 -D -S -gu -f ../reference/ce10/ce10.fa -r {} %s | bcftools view -bvcg - > ../vcf/raw.%s.{}.bcf'" % (' '.join(fasta_list),  outfile)

os.system(com)
# awk '{ if ($6>=%s ||  $1 ~ /^#/) print}' - |
os.system("cat ../ancillary/chr3ranges.txt | bcftools cat `ls -v ../vcf/raw.%s.*.bcf` | bcftools view - | vcf-sort > ../vcf/%s.vcf" % (outfile , outfile))
#os.system("bcftools index ../vcf/%s.bcf > ../vcf/%s.bcf" % (outfile, outfile))

# Output in multiple qualities
for quality in [10,20,30,40]:
	os.system("vcfutils.pl varFilter -d 3 -Q %s ../vcf/%s.vcf > ../vcf/%s.Q%s.vcf" % (quality, outfile, outfile, quality))
	# Zip, and index with tabix
	os.system("bgzip -f ../vcf/%s.Q%s.vcf" % (outfile, quality))
	os.system("tabix -f ../vcf/%s.Q%s.vcf.gz" % (outfile, quality))

# Remove temporary files
os.system("rm -f ../vcf/raw.%s.*" % outfile)



# Print total time
f.write("Time to complete: %s\n\n" % str(datetime.now()-startTime))