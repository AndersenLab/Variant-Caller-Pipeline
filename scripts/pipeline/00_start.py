#!/usr/bin/python

#=============#
# Description #
#=============#
#
# This script will submit a job that:
#
# * Aligns, creates bams
# * Sorts
# * Indexes

#SBATCH --job-name=indseq

#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=16384


#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/fq


#=================#
# Process Fastq's #
#=================#

# First - check whether the files in question exist.
import os, glob, sys, subprocess, re, gzip, math, tempfile

sys.path.append("/lscr2/andersenlab/dec211/scripts/pipeline/")
#sys.path.append("/exports/people/andersenlab/dec211/python_modules/peewee/")
#sys.path.remove("/usr/local/lib/python2.6/dist-packages/peewee-2.3.2-py2.6.egg")

from seq_utility_functions import *
from conn import *
from datetime import datetime
from peewee import *


#=======#
# Setup #
#=======#

# Determine system type:
system_type =  os.uname()[0]

# If local cd into appropriate directory
if system_type == "Darwin":
	os.chdir("../../data/fq")
	db_loc = "../"
else:
	db_loc = "/exports/people/andersenlab/dec211/"

# Setup Scripts Dir
scripts_dir = "../../scripts/pipeline/"

# Get Arguments
strain = sys.argv[1]

process_steps = {
		"debug_sqlite" : True,    # Uses an alternative database.
		"md5" : True,            # Runs an MD5 hash on every file and saves to database
		"fastq_stats" : True,     # Produce fastq stats on number of reads, unique reads, etc.
		"align" : True,
		"bam_stats" : True,
		"call_variants" : True,

}

#=========#
# Options #
#=========#

# Determine system type:
system_type =  os.uname()[0]
reference = "c_elegans.WS235.genomic"
md5_system = {"Linux" : "md5sum", "Darwin":"md5"} # Needed to make md5 work locally and on cluster.
system_cores = {"Linux": 8, "Darwin": 4} # For use with BWA

## LCR File
LCR_File = "WS235.wormbase.masked.bed.gz"
genome_length = int(file("../reference/%s/%s.fa.gz.amb" % (reference, reference), 'r').read().strip().split(" ")[0])

#=================#
# Process Fastq's #
#=================#

fq_set = glob.glob("*-%s-*fq.gz" % strain)

print fq_set

# Do we have fqs
if (len(fq_set) == 0):
	raise Exception("Error - No files with given name exist")

# Do we have pairs for each fq?
fq_pairs_count = ['-'.join(x.split("-")[0:3]) for x in fq_set]
if all([fq_pairs_count.count(x)==2 for x in fq_pairs_count]) != True:
	raise Exception("Pair Error")

# Fastq pairs
fastq_pairs = zip(sorted([x for x in fq_set if x.find("1.fq.gz") != -1]), sorted([x for x in fq_set if x.find("2.fq.gz") != -1]))

# Generate MD5's for each fastq
if process_steps["md5"] == True:
	save_md5(fq_set, strain)

#=================#
# Align Genome    #
#=================#

bam_set = [] # Created bams will later be merged.

# Assign sample using original fastq.
sample = fq_set[0].split("-")[2]

# Align genomes using bwa.
if process_steps["align"] == True:
	for fastq_set in fastq_pairs:
		# Generate bam name
		reference_loc = "../reference/%s/%s.fa.gz" % (reference, reference)
		bam_name = '-'.join(fastq_set[0].split("-")[0:4] + fastq_set[1].split("-")[3:4])
		split_bam_name = bam_name.split("-")
		library_LB = split_bam_name[1]
		readEntity_Group = '@RG\\tID:%s\\tLB:%s\\tSM:%s\\tPL:ILLUMINA' % (bam_name, library_LB, sample)

		# Align, add read-Entity_Group (@RG) header line, and output as a bam.
		os.system("bwa mem -M -t %s -R \"%s\" %s %s %s | samtools view -b -S -h -@ 2 -  > %s.tmp.bam" % (system_cores[system_type], readEntity_Group, reference_loc, fastq_set[0], fastq_set[1],bam_name))
		os.system("samtools sort -O bam -T sorting -@ %s %s.tmp.bam > %s.tmp.sorted.bam" % (system_cores[system_type], bam_name, bam_name))

		## Mark Duplicates, and remove.
		remove_duplicates = """
		mark_dups=`which MarkDuplicates.jar`
		java -jar $mark_dups \
		I=%s.tmp.sorted.bam \
		O=%s.bam \
		M=%s.duplicate_report.txt \
		VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true""" % (bam_name, bam_name, bam_name)
		os.system(remove_duplicates)

		# Parse Duplicate Stats From Picard and Store
		f = subprocess.check_output("egrep -v '^(#|^$)' %s.duplicate_report.txt" % bam_name, shell=True)
		dup_report = [x.split("\t") for x in f.split("\n")[0:2]]
		dup_report = zip([x.title() for x in dup_report[0]], dup_report[1])
		for record in dup_report:
			save_eav(strain, record[0], record[1], Entity_Group = "Duplication Statistics", Sub_Entity = bam_name + ".bam",  Tool="PICARD")

		# Index Bam
		os.system("samtools index %s.bam" % bam_name)


		
		# Store Bam Statistics
		samtools_stats = subprocess.check_output("samtools stats %s.bam | grep '^SN'" % bam_name , shell=True )
		for line in [x.split('\t')[1:3] for x in samtools_stats.split("\n")[:-1]]:
			save_eav(strain, line[0].replace(":","").title(), line[1],Entity_Group = "BAM Statistics", Sub_Entity = bam_name + ".bam", Tool = "Samtools Stats")

		# Generate Bam Depth and Coverage Statistics and save to database.
		for i in coverage(bam_name + ".bam"):
			print i, "ind"
			save_eav(strain, i[1], i[2], Sub_Attribute=i[0], Entity_Group = "BAM Statistics", Sub_Entity = bam_name + ".bam", Tool = "Samtools Depth + Python")
		
		# Generate md5 of bam and store
		save_md5([bam_name + ".bam"], Entity = strain)

		# Remove temporary files
		os.remove("%s.tmp.sorted.bam" % bam_name)
		os.remove("%s.tmp.bam" % bam_name)
		os.remove("%s.duplicate_report.txt" % bam_name)

		bam_set.append(bam_name + ".bam")

	# Merge BAMs into a single BAM.
	if len(bam_set) > 1:
		os.system("samtools merge -f -@ 4 %s.bam %s " % (strain, ' '.join(bam_set)))
	else:
		os.system("mv %s %s.bam" % (bam_set[0] , strain))
	# Create record indicating the bams that were merged.
	save_eav(strain, sample + ".bam", ', '.join(bam_set), Sub_Entity = "Constitutive Bams", Entity_Group = "Bam Merged Statistics")
	save_eav(strain, sample + ".bam", len(bam_set), Sub_Entity = "Constitutive Bam Count", Entity_Group = "Bam Merged Statistics")
	os.system("samtools index %s.bam" % sample)

# Generate Depth and coverage statistics of the merged bam
if process_steps["bam_stats"] == True:
	for i in coverage(sample + ".bam", "CHROMOSOME_MtDNA"):
		print i, "Merged"
		save_eav(strain, i[1], i[2], Sub_Entity = sample + ".bam", Sub_Attribute=i[0], Entity_Group = "Bam Merged Statistics", Tool = "Samtools Depth + Python")
		
	# Store Bam Statistics
	samtools_stats = subprocess.check_output( "samtools stats %s.bam | grep '^SN'" % sample , shell=True )
	line_dict = {}
	for line in [x.split('\t')[1:3] for x in samtools_stats.split("\n")[:-1]]:
		attribute = line[0].replace(":","").title()
		line_dict[attribute] = line[1] 
		save_eav(strain, attribute, line[1], Entity_Group = "BAM Merged Statistics", Sub_Entity = sample + ".bam", Tool = "Samtools")

	# Also Calculate Coverage as (# Reads * Avg. Read Length / Lenght Genome)
	calc_coverage = float(line_dict["Reads Mapped"]) * float(line_dict["Average Length"]) / genome_length
	save_eav(strain, "Depth of Coverage", calc_coverage, Entity_Group = "BAM Merged Statistics", Sub_Entity = sample + ".bam", Sub_Attribute = "(Reads Mapped * Avg. Read Length / Genome Length)", Tool = "Samtools")

# Remove intermediates.
for b in bam_set:
	try:
		pass
		#os.remove(b)
		#os.remove(b + ".bai")
	except:
		pass

# Finally - Generate and Save an md5 sum.
save_md5([sample + ".bam"], Entity = strain)

# Fetch contigs
header = subprocess.check_output("samtools view -H %s.bam" % strain, shell = True)
# Extract contigs from header and convert contigs to integers
contigs = {}
for x in re.findall("@SQ	SN:(?P<chrom>[A-Za-z0-9_]*)\WLN:(?P<length>[0-9]+)", header):
	contigs[x[0]] = int(x[1])

# Move bam to bam folder.
os.system("mv %s.bam ../bam/%s.bam" % (strain, sample))
os.system("mv %s.bam.bai ../bam/%s.bam.bai" % (strain, sample))
os.chdir("../bam/") # Change Dirs into bam directory.

#====================#
# Variant Calling    # - For initial calling
#====================#

# Generate and store bcf stats
def gen_bcf_stats(strain, bcf, description, Sub_Attribute):
	# Generate stats - Following Het Polarization, Overwrites previous stat file.
	tmp =  tempfile.NamedTemporaryFile().name
	os.system("bcftools stats --samples - %s > %s" % (bcf, tmp))

	# Get variant statistics.
	for i in import_stats("%s" % (tmp)):
		save_eav(strain, i[0], i[1], Entity_Group = "BCF Statistics", Sub_Entity = "%s.%s.bcf" % (strain, description), Sub_Attribute = Sub_Attribute, Tool="bcftools")


if process_steps["call_variants"] == True:
	command = "parallel --gnu 'samtools mpileup -t DP,DV,DP4,SP -g -f ../reference/%s/%s.fa.gz -r {} %s.bam | bcftools call --format-fields GQ,GP -c -v > ../individual_bcf/raw.%s.{}.bcf' ::: %s" % (reference, reference, strain, strain, ' '.join(contigs.keys()))
	os.system(command)
	os.system("echo %s | bcftools concat -O b `ls -v ../individual_bcf/raw.%s.*.bcf` > ../individual_bcf/%s.single.bcf" % (contigs.keys(), sample , sample.replace(".txt","")))
	
	# Index
	os.system("bcftools index -f ../individual_bcf/%s.single.bcf" % sample)

	# Remove temporary files
	os.system("rm -f ../individual_bcf/raw.%s.*" % sample)
