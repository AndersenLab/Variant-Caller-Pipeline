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

# Parse arguments

import sys, os

options = sys.argv

#==============#
# Depth Filter #
#==============#

if "-d" in options:
	depth_filter = options[options.index("-d") + 1]
else:
	depth_filter = ""
print depth_filter


# Output in multiple qualities
#for quality in [40]:
#	os.system("vcfutils.pl varFilter -d 3 -Q %s ../vcf/%s.vcf > ../vcf/%s.Q%s.vcf" % (quality, outfile, outfile, quality))
#	# Zip, and index with tabix
#	os.system("bgzip -f ../vcf/%s.Q%s.vcf" % (outfile, quality))
#	os.system("tabix -f ../vcf/%s.Q%s.vcf.gz" % (outfile, quality))