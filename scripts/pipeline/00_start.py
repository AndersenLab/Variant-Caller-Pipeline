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
from db_conn import *
from datetime import datetime
from peewee import *


def file_greater_than_0(filename):
	if os.path.isfile(filename) and os.path.getsize(filename) > 0:
		return True
	else:
		return False

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
for fastq_set in fastq_pairs:
	# Generate bam name
	reference_loc = "../reference/%s/%s.fa.gz" % (reference, reference)
	bam_name = '-'.join(fastq_set[0].split("-")[0:4] + fastq_set[1].split("-")[3:4])
	split_bam_name = bam_name.split("-")
	library_LB = split_bam_name[1]
	readEntity_Group = '@RG\\tID:%s\\tLB:%s\\tSM:%s\\tPL:ILLUMINA' % (bam_name, library_LB, sample)

	# Check that bam does not already exist
	if not file_greater_than_0("../bam_individual/" + bam_name + ".bam"):
		# Align, add read-Entity_Group (@RG) header line, and output as a bam.
		os.system("bwa mem -M -t %s -R \"%s\" %s %s %s | samtools view -b -S -h -@ 2 -  > ../bam_individual/%s.tmp.bam" % (system_cores[system_type], readEntity_Group, reference_loc, fastq_set[0], fastq_set[1],bam_name))
		os.system("samtools sort -O bam -T sorting -@ %s ../bam_individual/%s.tmp.bam > ../bam_individual/%s.tmp.sorted.bam" % (system_cores[system_type], bam_name, bam_name))

		## Mark Duplicates, and remove.
		remove_duplicates = """
		mark_dups=`which MarkDuplicates.jar`
		java -jar $mark_dups \
		I=../bam_individual/%s.tmp.sorted.bam \
		O=../bam_individual/%s.bam \
		M=../bam_individual/%s.duplicate_report.txt \
		VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true""" % (bam_name, bam_name, bam_name)
		os.system(remove_duplicates)

		# Parse Duplicate Stats From Picard and Store
		f = subprocess.check_output("egrep -v '^(#|^$)' ../bam_individual/%s.duplicate_report.txt" % bam_name, shell=True)
		dup_report = [x.split("\t") for x in f.split("\n")[0:2]]
		dup_report = zip([x.title() for x in dup_report[0]], dup_report[1])
		for record in dup_report:
			save_eav(strain, record[0], record[1], Entity_Group = "Duplication Statistics", Sub_Entity = bam_name + ".bam",  Tool="PICARD")

		# Index Bam
		os.system("samtools index ../bam_individual/%s.bam" % bam_name)


		
		# Store Bam Statistics
		samtools_stats = subprocess.check_output("samtools stats ../bam_individual/%s.bam | grep '^SN'" % bam_name , shell=True )
		for line in [x.split('\t')[1:3] for x in samtools_stats.split("\n")[:-1]]:
			save_eav(strain, line[0].replace(":","").title(), line[1],Entity_Group = "BAM Statistics", Sub_Entity = bam_name + ".bam", Tool = "Samtools Stats")

		# Generate Bam Depth and Coverage Statistics and save to database.
		for i in coverage("../bam_individual/" + bam_name + ".bam"):
			save_eav(strain, i[1], i[2], Sub_Attribute=i[0], Entity_Group = "BAM Statistics", Sub_Entity = bam_name + ".bam", Tool = "Samtools Depth + Python")
		
		# Generate md5 of bam and store
		save_md5(["../bam_individual/" + bam_name + ".bam"], Entity = strain)

		# Remove temporary files
		os.remove("../bam_individual/%s.tmp.sorted.bam" % bam_name)
		os.remove("../bam_individual/%s.tmp.bam" % bam_name)
		os.remove("../bam_individual/%s.duplicate_report.txt" % bam_name)

		bam_set.append("../bam_individual/" + bam_name + ".bam")
	else:
		bam_set = glob.glob("../bam_individual/*-%s-*bam" % strain)



# Merge BAMs into a single BAM.
if not file_greater_than_0("../bam_merged/%s.bam" % (strain)):
	if len(bam_set) > 1:
		os.system("samtools merge -f -@ 6 ../bam_merged/%s.bam %s" % (strain, ' '.join(bam_set)))
	else:
		os.system("mv ../bam_individual/%s ../bam_merged/%s.bam" % (bam_set[0], strain))
	# Create record indicating the bams that were merged.
	save_eav(strain, sample + ".bam", ', '.join(bam_set), Sub_Entity = "Constitutive Bams", Entity_Group = "Bam Merged Statistics")
	save_eav(strain, sample + ".bam", len(bam_set), Sub_Entity = "Constitutive Bam Count", Entity_Group = "Bam Merged Statistics")

# Index as needed
if not file_greater_than_0("../bam_merged/%s.bam.bai" % (strain)):
	os.system("samtools index ../bam_merged/%s.bam" % sample)

# Generate Depth and coverage statistics of the merged bam
# Test to see if they have already been saved; 77 bam stats are saved/merged bam
if EAV.select().where(EAV.Sub_Entity == strain).count() != 77:
	for i in coverage("../bam_merged/" + sample + ".bam", "CHROMOSOME_MtDNA"):
		save_eav(strain, i[1], i[2], Sub_Entity = sample + ".bam", Sub_Attribute=i[0], Entity_Group = "Bam Merged Statistics", Tool = "Samtools Depth + Python")
		
	# Store Bam Statistics
	samtools_stats = subprocess.check_output( "samtools stats ../bam_merged/%s.bam | grep '^SN'" % sample , shell=True )
	line_dict = {}
	for line in [x.split('\t')[1:3] for x in samtools_stats.split("\n")[:-1]]:
		attribute = line[0].replace(":","").title()
		line_dict[attribute] = line[1] 
		save_eav(strain, attribute, line[1], Entity_Group = "Bam Merged Statistics", Sub_Entity = sample + ".bam", Tool = "Samtools")

	# Also Calculate Coverage as (# Reads * Avg. Read Length / Lenght Genome)
	calc_coverage = float(line_dict["Reads Mapped"]) * float(line_dict["Average Length"]) / genome_length
	save_eav(strain, "Depth of Coverage", calc_coverage, Entity_Group = "BAM Merged Statistics", Sub_Entity = sample + ".bam", Sub_Attribute = "(Reads Mapped * Avg. Read Length / Genome Length)", Tool = "Samtools")


	# Finally - Generate and Save an md5 sum.
	save_md5(["../bam_merged/" + sample + ".bam"], Entity = strain)

# Fetch contigs
header = subprocess.check_output("samtools view -H ../bam_merged/%s.bam" % strain, shell = True)
# Extract contigs from header and convert contigs to integers
contigs = {}
for x in re.findall("@SQ	SN:(?P<chrom>[A-Za-z0-9_]*)\WLN:(?P<length>[0-9]+)", header):
	contigs[x[0]] = int(x[1])


#====================#
# Variant Calling    # - For initial (individual) calling
#====================#

# Generate and store bcf stats
def gen_bcf_stats(strain, bcf, description, Sub_Attribute):
	# Generate stats - Following Het Polarization, Overwrites previous stat file.
	tmp =  tempfile.NamedTemporaryFile().name
	os.system("bcftools stats --samples - %s > %s" % (bcf, tmp))

	# Get variant statistics.
	for i in import_stats("%s" % (tmp)):
		save_eav(strain, i[0], i[1], Entity_Group = "BCF Statistics", Sub_Entity = "%s.%s.bcf" % (strain, description), Sub_Attribute = Sub_Attribute, Tool="bcftools")

if not file_greater_than_0("../individual_bcf/%s.single.bcf" % sample):
	command = "echo %s | xargs -P 6 -I {} sh -c'samtools mpileup -t DP,DV,DP4,SP -g -f ../reference/%s/%s.fa.gz -r {} ../bam_merged/%s.bam | bcftools call --format-fields GQ,GP -c -v > ../individual_bcf/raw.%s.{}.bcf' ::: %s" % (' '.join(contigs.keys(), reference, reference, strain, strain))
	os.system(command)
	os.system("echo %s | bcftools concat -O b `ls -v ../individual_bcf/raw.%s.*.bcf` > ../individual_bcf/%s.single.bcf" % (contigs.keys(), sample , sample.replace(".txt","")))

	# Index
	os.system("bcftools index -f ../individual_bcf/%s.single.bcf" % sample)

	# Remove temporary files
	os.system("rm -f ../individual_bcf/raw.%s.*" % sample)
