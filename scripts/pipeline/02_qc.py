#!/usr/bin/python

#=============#
# Description #
#=============#
#
# This script will submit a job that:
#
# * Runs reports
#



#SBATCH --job-name=vcfcomp

#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem=16384

#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/vcf

#---------#
# Imports #
#---------#

import os, sys
import glob
from itertools import combinations
from datetime import datetime

print "Job Id: %s\n" % os.system("echo $SLURM_JOB_ID")
print "Command: %s\n" % (sys.argv)

# Determine type of qc being performed. 1 file, 2 files (compare) or comparison among many files (multiple).

# For 1 file, concordance and depth stats are produced for all pairwise possible.

for vcf in combinations(glob.glob('*%s*.vcf.gz' % sys.argv[1]),2):
		outfile = '_'.join([vcf[0].replace(".vcf.gz",""),vcf[1].replace(".vcf.gz","")])
		os.system('mkdir ../reports/%s' % sys.argv[1])
		os.system('bcftools gtcheck -s %s %s > ../reports/%s/%s.txt' % (vcf[0],vcf[1], sys.argv[1], outfile))
		os.system('python ../../scripts/analysis/custom_gtcheck.py ../reports/%s/%s_%s.txt' % (outfile,vcf[0],vcf[1]))


"""
os.system('bcftools stats -s - %s > ../reports/%s.stats.txt' % (" ".join(sys.argv[2:]), outfile))

# Genotype check
p=file("../reports/plot_commands.txt","a+")
p.write("scp dec211@dhunni.biochem.northwestern.edu:/lscr2/andersenlab/dec211/data/reports/%s.stats.txt %s.stats.txt\n" % (outfile, outfile))
p.write("scp dec211@dhunni.biochem.northwestern.edu:/lscr2/andersenlab/dec211/data/vcf/%s %s\n" % (sys.argv[2],sys.argv[2]))
p.write("scp dec211@dhunni.biochem.northwestern.edu:/lscr2/andersenlab/dec211/data/vcf/%s %s\n" % (sys.argv[3],sys.argv[3]))
p.write("scp dec211@dhunni.biochem.northwestern.edu:/lscr2/andersenlab/dec211/data/vcf/%s.tbi %s.tbi\n" % (sys.argv[2],sys.argv[2]))
p.write("scp dec211@dhunni.biochem.northwestern.edu:/lscr2/andersenlab/dec211/data/vcf/%s.tbi %s.tbi\n" % (sys.argv[3],sys.argv[3]))
p.write('bcftools gtcheck -H -p %s/ -g %s \n' % (outfile, " ".join(sys.argv[2:])))
p.write("plot-vcfstats --main-title %s -t %s -t %s -r -p %s/ %s.stats.txt\n" % (outfile, sys.argv[2].replace(".txt","").replace(".vcf.gz",""), sys.argv[3].replace(".txt","").replace(".vcf.gz",""),outfile, outfile))
p.write("\n\n")
"""
# Print total time

