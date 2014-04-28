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
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=16384

#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/vcf



from datetime import datetime
startTime = datetime.now()
import os, sys

outfile = sys.argv[1].replace('[','\[').replace(']','\]')
infiles = [x.replace('[','\[').replace(']','\]') for x in sys.argv[2:]]

job_id = os.system("echo $SLURM_JOB_ID")


print "Job Id: %s" % job_id
print
print "Command: %s" % (sys.argv)
print

os.system('bcftools stats %s > ../reports/%s.txt' % (outfile, infiles))
os.system('plot-vcfstats ../reports/%s.txt -p %s/' % (outfile, outfile))

# Print total time
print("Time to complete: %s" % str(datetime.now()-startTime))
