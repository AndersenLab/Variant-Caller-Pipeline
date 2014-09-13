#!/bin/python 

#=============#
# Description #
#=============#
#
# This script will submit a job that:
#
# * Aligns, creates bams
# * Sorts
# * Indexes

#SBATCH --job-name=bwa

#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=16384


#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data


#=================#
# Process Fastq's #
#=================#

# First - check whether the files in question exist.
import os, glob, sys, subprocess
import sys
import re
import gzip
from seq_utility_functions import *
from datetime import datetime
from peewee import *

#=======#
# Setup #
#=======#

# Determine system type:
system_type =  os.uname()[0]

# The dictionary process_steps is used to skip steps for debugging purposes.
process_steps = {
	"debug_sqlite" : True,    # Uses an alternative database.
	"md5" : False,            # Runs an MD5 hash on every file and saves to database
	"fastq_stats" : True,     # Produce fastq stats on number of reads, unique reads, etc.
}

#==========#
# Database #
#==========#

# Drop tables if specified.
if process_steps["debug_sqlite"] == True:
	db = SqliteDatabase('DEBUG_seq_data.db')
else:
	db = SqliteDatabase("seq_data.db").connect()

db.connect()


class EAV(Model):
    Group = CharField(index=True)
    Entity = CharField(index=True)
    Sub_Entity = CharField(index=True)
    Attribute = CharField(index=True)
    Value = CharField(index=True)
    Tool = CharField(index=True, null=True)
    DateTime = DateTimeField()

    class Meta:
        database = db # This model uses the "seq_data.db" database.
        indexes = (
            # Specify a unique multi-column index on from/to-user.
            (('Group', 'Entity', 'Attribute', 'Value'), True),
        )



# Drop tables if specified.
if process_steps["debug_sqlite"] == True:
	db.drop_tables([EAV], safe=True)
	db.create_tables([EAV], safe=True)
else:
	db.create_tables([EAV], safe=True)


def store_eav(Group, Entity, Attribute, Value, Sub_Entity=None, Tool=None):
	record = EAV(Group=Group, Entity=Entity, Attribute=Attribute, Value=Value, DateTime=datetime.now(), Sub_Entity=Sub_Entity, Tool=Tool)
	try:
		record.save()
	except:
		pass

#===========#
# Functions #
#===========#

def calc_bam_depth_coverage(Group, reference, bam, by_chr = False):
	""" 
		Calculates the average Depth, Covered Bases, and Coverage of a given Bam File
	"""
	results = subprocess.check_output("samtools depth %s.bam | awk '{sum+=$3;cnt++}END{print sum/cnt \"\t\" sum}'" % bam.replace(".bam",""), shell=True).replace("\n", "").split("\t")
	store_eav(Group, bam + ".bam", "Average Depth", results[0], Tool = "Samtools + Awk")
	store_eav(Group, bam + ".bam", "Covered Bases", results[1], Tool = "Samtools + Awk")

	# Calculate Coverage - ALL and by Chromosome
	reference_length = file("../reference/%s/%s.fa.amb" % (reference, reference), 'r').read().split(" ")[0]
	store_eav(Group, bam + ".bam", "Coverage", int(results[1])/float(reference_length), Sub_Entity = "Genome", Tool = "Samtools + Awk")

	# Calculate Average Depth and Coverage for each chromosome indivudally.
	if by_chr == True:
		# For merged bam, generate per-chromosome stats.
		contigs = [x.split("\t")[:2] for x in file("../reference/%s/%s.fa.fai" % (reference, reference), 'r').read().split("\n")[:-1]]
		for chr in contigs:
			results = subprocess.check_output("samtools depth -r %s %s.bam | awk '{sum+=$3;cnt++}END{print sum/cnt \"\t\" sum}'" % (chr[0], bam.replace(".bam","")), shell=True).replace("\n", "").split("\t")
			store_eav(Group, bam + ".bam", "Average Depth", results[0], Sub_Entity = chr[0],Tool = "Samtools + Awk")
			store_eav(Group, bam + ".bam", "Covered Bases", results[1], Sub_Entity = chr[0], Tool = "Samtools + Awk")


#=========#
# Options #
#=========#

reference = "ce10"
md5_system = {"Linux" : "md5sum", "Darwin":"md5"} # Needed to make md5 work locally and on cluster.
system_cores = {"Linux": 8, "Darwin": 4} # For use with BWA

#=================#
# Process Fastq's #
#=================#

args = sys.argv[1:] # Used to specify a strain name for fq's that will be aligned.

fq_set = glob.glob("*-%s-*fq.gz" %args[0])

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
	md5 = subprocess.check_output("parallel %s ::: %s" % (md5_system[system_type], ' '.join(fq_set)), shell=True)
	m = re.findall('MD5 \((.*)\) = (.*)', md5)
	for i in m:
		# Check File Hashes
		if i[0].split("-")[3] != i[1][0:5]:
			raise Exception("MD5 Hash for %s Mismatch. Hash for file is %s" % (i[0],i[1]))
		# Store if Correct
		store_eav("Fastq Statistics", i[0], "File Hash", i[1], Tool="MD5")


# Generate sequence fastq statistics using awk
command = "parallel \"gunzip -c {} | awk -v filename={} '((NR-2)%%4==0){read=\$1; total++; count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print filename,total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'\" ::: %s" % ' '.join(fq_set)

if process_steps["fastq_stats"] == True:
	fq_stats = subprocess.check_output(command, shell=True)
	output = subprocess.check_output(command, shell=True)

	for i in [x.split(" ") for x in fq_stats.split("\n")[:-1]]:
		# Save fastq statistics
		# i[0] -> Filename
		store_eav("Fastq Statistics", i[0], "Number of Reads", i[1], Tool="Awk")
		store_eav("Fastq Statistics", i[0], "Unique Reads", i[2], Tool="Awk")
		store_eav("Fastq Statistics", i[0], "Frequency of Unique Reads", i[3], Tool="Awk")
		store_eav("Fastq Statistics", i[0], "Most abundant Sequence", i[4], Tool="Awk")
		store_eav("Fastq Statistics", i[0], "Number of times most abundant sequence occurs", i[5], Tool="Awk")
		store_eav("Fastq Statistics", i[0], "Frequency of Most Abundant Sequence", i[6], Tool="Awk")

for fq in fq_set:
	fastq_info = extract_fastq_info(fq)
	store_eav("Fastq Meta", fq, "Flowcell Lane", fastq_info["flowcell_lane"], Tool="Python Script")
	store_eav("Fastq Meta", fq, "Index", fastq_info["index"], Tool="Python Script")
	store_eav("Fastq Meta", fq, "Instrument", fastq_info["instrument"][1:], Tool="Python Script")
	store_eav("Fastq Meta", fq, "Read Pair", fastq_info["pair"], Tool="Python Script")

#=================#
# Align Genome    #
#=================#

bam_set = [] # Created bams will later be merged.

# Align genomes using bwa.
for fastq_set in fastq_pairs:

	# Generate bam name
	reference_loc = "../reference/%s/%s.fa" % (reference, reference)
	bam_name = '-'.join(fastq_set[0].split("-")[0:4] + fastq_set[1].split("-")[3:4])
	split_bam_name = bam_name.split("-")
	library_LB = split_bam_name[1]
	sample = split_bam_name[2]
	readgroup = '@RG\\tID:%s\\tLB:%s\\tSM:%s\\tPL:ILLUMINA' % (bam_name, library_LB, sample)

	# Align, add read-group (@RG) header line, and output as a bam.
	os.system("bwa mem -M -t %s -R \"%s\" %s %s %s | samtools view -b -S -h -@ 2 -  > %s.tmp.bam" % (system_cores[system_type], readgroup, reference_loc, fastq_set[0], fastq_set[1],bam_name))
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
	print "egrep -v '^(#|^$)' %s.duplicate_report.txt" % bam_name
	f = subprocess.check_output("egrep -v '^(#|^$)' %s.duplicate_report.txt" % bam_name, shell=True)
	dup_report = [x.split("\t") for x in f.split("\n")[0:2]]
	dup_report = zip([x.title() for x in dup_report[0]], dup_report[1])

	for record in dup_report:
		store_eav("BAM Statistics", bam_name + ".bam", record[0], record[1], Tool="PICARD")

	# Store Bam Statistics
	samtools_stats = subprocess.check_output( "samtools stats %s.bam | grep '^SN'" % bam_name , shell=True )
	for line in [x.split('\t')[1:3] for x in samtools_stats.split("\n")[:-1]]:
		store_eav("BAM Statistics", bam_name + ".bam", line[0].replace(":","").title(), line[1], Tool = "Samtools")

	# Generate Bam Depth and Coverage Statistics and save to database.
	calc_bam_depth_coverage("BAM Statistics", reference, bam_name)

	# Generate md5 of bam and store
	md5 = subprocess.check_output("%s %s.bam" % (md5_system[system_type], bam_name), shell=True)
	m = re.match('MD5 \((.*)\) = (.*)', md5)
	store_eav("BAM Statistics", m.group(1), "File Hash", m.group(2), Tool="MD5")

	# Remove temporary files
	os.remove("%s.tmp.sorted.bam" % bam_name)
	os.remove("%s.tmp.bam" % bam_name)
	os.remove("%s.duplicate_report.txt" % bam_name)

	bam_set.append(bam_name + ".bam")

# Merge BAMs into a single BAM.
os.system("samtools merge -f -@ 4 %s.bam %s " % (sample, ' '.join(bam_set)))
os.system("samtools index %s.bam" % sample)

# Generate Depth and coverage statistics of the merged bam
calc_bam_depth_coverage("BAM Merged Statistics", reference, sample, by_chr = True)
# Store Bam Statistics
samtools_stats = subprocess.check_output( "samtools stats %s.bam | grep '^SN'" % sample , shell=True )
for line in [x.split('\t')[1:3] for x in samtools_stats.split("\n")[:-1]]:
	store_eav("BAM Merged Statistics", sample + ".bam", line[0].replace(":","").title(), line[1], Tool = "Samtools")

# Index Merged Bam File, and remove intermediates.
for b in bam_set:
	os.remove(b)

# Finally - Generate and Save an md5 sum.
md5 = subprocess.check_output("%s %s.bam" % (md5_system[system_type], sample), shell=True)
m = re.match('MD5 \((.*)\) = (.*)', md5)
store_eav("BAM Merged Statistics", m.group(1), "File Hash", m.group(2), Tool="MD5")

#======================================#
# Generate Mitochondrial Statistics    #
#======================================#

