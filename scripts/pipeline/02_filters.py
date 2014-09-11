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

import sys, os

vcf = sys.argv[1]
options = sys.argv[2:]

print "VCF", vcf
print "Options", options
#==============#
# Depth Filter #
#==============#

filter_set = []

# Depth
if "-d" in options:
	d_option = options[options.index("-d") + 1]
	if d_option == "avg2":
		threshold = int(float(os.popen("bcftools query -f '%%DP\\n' %s | awk '{ total += $1; count++ } END { avg=(total/count);  print (avg + 3*sqrt(avg)) }'" % vcf).read().strip()))
	elif d_option == "avg3":
		threshold = int(float(os.popen("bcftools query -f '%%DP\\n' %s | awk '{ total += $1; count++ } END { avg=(total/count);  print (avg + 2*sqrt(avg)) }'" % vcf).read().strip()))
	else:
		threshold = int(d_option)
	# Set up depth filter
	filter_set.append("bcftools filter -O b --include 'DP<%s' " % (threshold))
	d_file = ".d%04d" % (threshold)
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
	filter_set.append("bcftools filter -O b --include '%%QUAL>%s' " % (q_option))
	q_file = ".Q%04d" % int(q_option)
else:
	q_file = ""

# Low Complexity Regions
if "-l" in options:
	lcr_filter = "-R ../ancillary/ce10_no_LCR.bed -O b "
	lcr_file = ".lcr"
else:
	lcr_filter = ""
	lcr_file = ""

# Multiallelic filter
if "-m" in options:
	filter_set.append("bcftools view -m2 -M2 -O b")
	m_file = ".m"
else:
	m_file = ""

print filter_set, options

# Set filter_variable
if len(filter_set) > 0:
	filter_set = " | " +  " | ".join(filter_set)
else:
	filter_set = ""


vcf_filename = ''.join([vcf.replace(".vcf.gz","").replace(".bcf","").replace(".vcf",""),".filter",q_file,d_file,lcr_file,m_file,".bcf"])



print "bcftools view %s %s %s > %s " % (lcr_filter, vcf, filter_set, vcf_filename)

os.system("bcftools view %s %s %s > %s " % (lcr_filter, vcf, filter_set, vcf_filename))
os.system("bcftools index -f %s" % vcf_filename)

