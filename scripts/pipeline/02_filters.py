#!/usr/bin/python

#=============#
# Description #
#=============#
#
# This script will submit a job that:
#
# * Filters VCF Files
#



#SBATCH --job-name=vcfFilt

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

import sys, os, subprocess

vcf = sys.argv[1]
options = sys.argv[2:]

#==============#
# Depth Filter #
#==============#

filter_set = []

# Depth
if "-d" in options:
	d_option = options[options.index("-d") + 1]
	if d_option == "avg2":
		threshold = float(subprocess.check_output("bcftools query -f '%%DP\n' %s | awk '{ total += $1; count++ } END { avg=(total/count);  print (avg + 3*sqrt(avg)) }'" % (vcf), shell=True))
	elif d_option == "avg3":
		threshold = float(subprocess.check_output("bcftools query -f '%%DP\n' %s | awk '{ total += $1; count++ } END { avg=(total/count);  print (avg + 2*sqrt(avg)) }'" % (vcf), shell=True))
	else:
		threshold = d_option
	# Set up depth filter
		filter_set.append("bcftools filter --include 'DP<%s' --soft-filter 'DP_lt_%s'" % (threshold , threshold))
		d_file = ".d%s" % threshold
else:
	d_file = ""

# Polarize Hets
if "-h" in options:
	polarize_hets = ""
	h_file = ".h"
else:
	polarize_hets = ""
	h_file = ""

# Quality
if "-Q" in options:
	q_option = options[options.index("-Q") + 1]
	filter_set.append("bcftools filter -O b --include '%%QUAL>%s' --soft-filter 'QUAL_gt_%s'" % (q_option, q_option))
	q_file = ".Q%s" % q_option
else:
	q_file = ""

# Low Complexity Regions
if "-l" in options:
	lcr_filter = "-R ../ancillary/ce10_no_LCR.bed"
	lcr_file = ".lcr"
else:
	lcr_filter = ""
	lcr_file = ""

vcf_filename = ''.join([vcf.replace(".vcf.gz","").replace(".bcf","").replace(".vcf",""),q_file,d_file,lcr_file,".bcf"])


# Set filter_variable
if "-d" in options or "-Q" in options:
	filter_set = " | ".join(filter_set)
else:
	filter_set = ""

print "bcftools view %s %s | %s > %s " % (lcr_filter, vcf, filter_set, vcf_filename)

os.system("bcftools view %s %s | %s > %s " % (lcr_filter, vcf, filter_set, vcf_filename))
os.system("bcftools index -f %s" % vcf_filename)

