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

options = sys.argv[1:]


#==============#
# Depth Filter #
#==============#

# Depth
if "-d" in options:
	depth_filter = options[options.index("-d") + 1]
	d_file = "d%s" % depth_filter
else:
	depth_filter = ""
	d_file = ""

# Polarize Hets
if "-h" in options:
	polarize_hets = ""
	h_file = ".h"
else:
	polarize_hets = ""
	h_file = ""

# Low Complexity Regions
if "-l" in options:
	lcr_filter = ""
	lcr_file = ".lcr"
else:
	lcr_filter = ""
	lcr_file = ""

os.system("vcfutils.pl varFilter -d 3 -D %s -Q 30 %s > %s.%s.vcf.gz" % (depth_filter, options[0], options[0].replace(".txt.vcf.gz", ""), d_file))

# Output in multiple qualities
#for quality in [40]:
#	os.system("vcfutils.pl varFilter -d 3 -Q %s ../vcf/%s.vcf > ../vcf/%s.Q%s.vcf" % (quality, outfile, outfile, quality))
#	# Zip, and index with tabix
#	os.system("bgzip -f ../vcf/%s.Q%s.vcf" % (outfile, quality))
#	os.system("tabix -f ../vcf/%s.Q%s.vcf.gz" % (outfile, quality))