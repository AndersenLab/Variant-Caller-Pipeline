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
system_type = subprocess.check_output("uname", shell=True).replace("\n","")

# The dictionary process_steps is used to skip steps for debugging purposes.
process_steps = {
	"debug_sqlite" : True,    # Uses an alternative database.
	"md5" : False,            # Runs an MD5 hash on every file and saves to database
	"fastq_stats" : False,     # Produce fastq stats on number of reads, unique reads, etc.
}

#=========#
# Options #
#=========#

reference = "ce10"
md5_system = {"Linux" : "md5sum", "Darwin":"md5"} # Needed to make md5 work locally and on cluster.
system_cores = {"Linux": 8, "Darwin": 4} # For use with BWA

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
    Attribute = CharField(index=True)
    Value = CharField(index=True)
    DateTime = DateTimeField()

    class Meta:
        database = db # This model uses the "seq_data.db" database.
        indexes = (
            # Specify a unique multi-column index on from/to-user.
            (('Group', 'Entity', 'Attribute', 'Value'), True),
        )

# Drop tables if specified.
if process_steps["debug_sqlite"] == True:
	#db.drop_tables([EAV], safe=True)
	db.create_tables([EAV], safe=True)
else:
	db.create_tables([EAV], safe=True)


def store_eav(Group, Entity, Attribute, Value):
	record = EAV(Group=Group, Entity=Entity, Attribute=Attribute, Value=Value, DateTime=datetime.now())
	try:
		record.save()
	except:
		pass

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
		store_eav("Fastq Statistics", i[0], "File Hash", i[1])
		if i[0].split("-")[3] != i[1][0:5]:
			raise Exception("MD5 Hash for %s Mismatch. Hash for file is %s" % (i[0],i[1]))


# Generate sequence fastq statistics using awk
command = "parallel \"gunzip -c {} | awk -v filename={} '((NR-2)%%4==0){read=\$1; total++; count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print filename,total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'\" ::: %s" % ' '.join(fq_set)

if process_steps["fastq_stats"] == True:
	fq_stats = subprocess.check_output(command, shell=True)
	output = subprocess.check_output(command, shell=True)

	for i in [x.split(" ") for x in fq_stats.split("\n")[:-1]]:
		# Save fastq statistics
		# i[0] -> Filename
		store_eav("Fastq Statistics", i[0], "Number of Reads", i[1])
		store_eav("Fastq Statistics", i[0], "Unique Reads", i[2])
		store_eav("Fastq Statistics", i[0], "Frequency of Unique Reads", i[3])
		store_eav("Fastq Statistics", i[0], "Most abundant Sequence", i[4])
		store_eav("Fastq Statistics", i[0], "Number of times most abundant sequence occurs", i[5])
		store_eav("Fastq Statistics", i[0], "Frequency of Most Abundant Sequence", i[6])

# Fetch data from fastq file (Machine, Barcode, Flowcell Lane, etc.)
fastq_info = extract_fastq_info("BGI1-RET6-JU1395-1f256-2.fq.gz")

for fq in fq_set:
	fastq_info = extract_fastq_info(fq)
	store_eav("Fastq Header", fq, "Flowcell Lane", fastq_info["flowcell_lane"])
	store_eav("Fastq Header", fq, "Index", fastq_info["index"])
	store_eav("Fastq Header", fq, "Instrument", fastq_info["instrument"])
	store_eav("Fastq Header", fq, "Pair", fastq_info["pair"])

#=================#
# Align Genome    #
#=================#

# Align genomes using bwa.
for fastq_set in fastq_pairs:
	# Generate bam name
	reference_loc = "../reference/%s/%s.fa" % (reference, reference)
	bam_name = '-'.join(fastq_set[0].split("-")[0:4] + fastq_set[1].split("-")[3:4])
	split_bam_name = bam_name.split("-")
	library_LB = split_bam_name[1]
	sample_SM = split_bam_name[2]
	readgroup = '@RG\\tID:%s\\tLB:%s\\tSM:%s\\tPL:ILLUMINA' % (bam_name, library_LB, sample_SM)

	# Align, add read-group (@RG) header line, and output as a bam.
	print("bwa mem -t %s -R \"%s\" %s %s %s | samtools view -b -S -h -@ 2 -  > %s.tmp.bam" % (system_cores[system_type], readgroup, reference_loc, fastq_set[0], fastq_set[1],bam_name))
	print "samtools sort -O bam -T sorting -@ %s > %s.tmp.sorted.bam" % (system_cores[system_type], bam_name)
	# Save Summary Stats
	#os.system("samtools stats %s | grep '^SN'" % bam_name)
	#samtools_stats = subprocess.check_output( samtools_stats, shell=True )

	## Mark Duplicates
	remove_duplicates = """java -jar MarkDuplicates.jar \
	I=%s.tmp.sorted.bam \
	O=%s.bam \
	M=${bam_name}.duplicate_report.txt \
	VALIDATION_STRINGENCY=SILENT \
	REMOVE_DUPLICATES=true"""


	# Store Bam Statistics
	#for line in [x.split('\t')[1:3] for x in samtools_stats.split("\n")[:-1]]:
	#	store_eav("BAM Statistics", bam_name, line[0].replace(":","").title(), line[1])

	# Generate md5 of bam and store
	#md5 = subprocess.check_output("%s %s.bam" % (md5_system[system_type], bam_name)
	#m = re.match('MD5 \((.*)\) = (.*)', md5)

	samtools_stats = """SN	raw total sequences:	1741156
SN	filtered sequences:	0
SN	sequences:	1741156
SN	is sorted:	0
SN	1st fragments:	870580
SN	last fragments:	870576
SN	reads mapped:	1733254
SN	reads mapped and paired:	1731612	# paired-end technology bit set + both mates mapped
SN	reads unmapped:	7902
SN	reads properly paired:	1726280	# proper-pair bit set
SN	reads paired:	1741156	# paired-end technology bit set
SN	reads duplicated:	0	# PCR or optical duplicate bit set
SN	reads MQ0:	104813	# mapped and MQ=0
SN	reads QC failed:	0
SN	non-primary alignments:	0
SN	total length:	156660257	# ignores clipping
SN	bases mapped:	155949077	# ignores clipping
SN	bases mapped (cigar):	154161626	# more accurate
SN	bases trimmed:	0
SN	bases duplicated:	0
SN	mismatches:	983937	# from NM fields
SN	error rate:	6.382503e-03	# mismatches / bases mapped (cigar)
SN	average length:	89
SN	maximum length:	90
SN	average quality:	64.3
SN	insert size average:	220.8
SN	insert size standard deviation:	54.0
SN	inward oriented pairs:	826481
SN	outward oriented pairs:	11444
SN	pairs with other orientation:	459
SN	pairs on different chromosomes:	1088
"""



"""

###############
# Get FASTQ's #
###############

# Convert f1 and f2 to arrays
f1=(${1//-/ })
f2=(${2//-/ })

echo $1
echo $2

# Set reference
ref=${3##*/}

bam_name=${BAM_DIR}${f1[0]}-${f1[1]}-${f1[2]}-${ref%.fa}-${f1[3]}-${f2[3]}.bam

bwa mem -t 8 ../reference/${3} ${1} ${2} | samtools view -b -S -h -  > ${BAM_DIR}${f1[3]}-${f2[3]}.tmp.bam

## Replace Header & Sort Simultaneously!
java	-jar	${PICARD}AddOrReplaceReadGroups.jar	\
I=${BAM_DIR}${f1[3]}-${f2[3]}.tmp.bam	\
O=${BAM_DIR}${f1[3]}-${f2[3]}.tmp.sorted.bam	\
RGID=`basename ${bam_name}`	RGLB=${f1[1]}	\
RGPL=ILLUMINA RGPU=${f1[1]} \
RGSM=${f1[2]}	RGCN=BGI \
SO=coordinate

# Remove temp
rm ${BAM_DIR}${f1[3]}-${f2[3]}.tmp.bam ${BAM_DIR}${f1[3]}-${f2[3]}.tmp.sorted.bam


# Index Bams
samtools index $bam_name
"""
