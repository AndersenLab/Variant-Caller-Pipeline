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

print sys.path

from seq_utility_functions import *
from datetime import datetime
from peewee import *


## Fix for python 2.6

if "check_output" not in dir( subprocess ): # duck punch it in!
    def f(*popenargs, **kwargs):
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd)
        return output
    subprocess.check_output = f

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
		"debug_sqlite" : False,    # Uses an alternative database.
		"md5" : True,            # Runs an MD5 hash on every file and saves to database
		"fastq_stats" : True,     # Produce fastq stats on number of reads, unique reads, etc.
		"align" : True,
		"bam_stats" : True,
		"call_variants" : True,

}

#=========#
# Options #
#=========#

reference = "ce10"
md5_system = {"Linux" : "md5sum", "Darwin":"md5"} # Needed to make md5 work locally and on cluster.
system_cores = {"Linux": 8, "Darwin": 4} # For use with BWA

## LCR File
LCR_File = "WS220.wormbase.masked.bed.gz"
genome_length = int(file("../reference/%s/%s.fa.amb" % (reference, reference), 'r').read().strip().split(" ")[0])

#==========#
# Database #
#==========#

# Drop tables if specified.
if process_steps["debug_sqlite"] == True:
	db = SqliteDatabase(db_loc + "DEBUG_seq_data.db")
else:
	db = SqliteDatabase(db_loc + "seq_data.db")

db.connect()

class EAV(Model):
    Entity_Group = CharField(index=True, null = True)
    Entity = CharField(index=True)
    Sub_Entity = CharField(index=True, null = True) # Used for option secondary super Entity_Groupings
    Attribute = CharField(index=True)
    Sub_Attribute = CharField(index=True, null = True) # Used for optional secondary sub Entity_Groupings
    Value = CharField(index=True)
    Tool = CharField(index=True, null=True)
    DateTime = DateTimeField()

    class Meta:
        database = db # This model uses the "seq_data.db" database.
        indexes = (
            # Specify a unique multi-column index on from/to-user.
            (('Entity', 'Attribute', 'Value',), True),
        )

# Drop tables if specified.
if process_steps["debug_sqlite"] == True:
	#db.drop_tables([EAV], safe=True)
	db.create_tables([EAV], safe=True)
else:
	db.create_tables([EAV], safe=True)


def save_eav(Entity, Attribute, Value, Entity_Group=None, Sub_Entity=None, Sub_Attribute=None, Tool=None):
	record = EAV(Entity=Entity, Attribute=Attribute, Value=Value, Entity_Group=Entity_Group, Sub_Entity=Sub_Entity, Sub_Attribute=Sub_Attribute, Tool=Tool, DateTime=datetime.now())
	try:
		record.save()
	except:
		record.update()

def save_md5(files = [], Entity = ""):
	md5 = subprocess.check_output("parallel %s ::: %s" % (md5_system[system_type], ' '.join(files)), shell=True)
	m = re.findall('MD5 \((.*)\) = (.*)', md5)
	for i in m:
		# Check File Hashes
		save_eav(Entity, "File Hash", i[1], Entity_Group = "MD5 Hash", Sub_Entity= i[0], Tool = "MD5")

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

# Generate sequence fastq statistics using awk
if process_steps["fastq_stats"] == True:
	command = "parallel \"gunzip -c {} | awk -v filename={} '((NR-2)%%4==0){read=\$1; total++; count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print filename,total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'\" ::: %s" % ' '.join(fq_set)

	if process_steps["fastq_stats"] == True:
		fq_stats = subprocess.check_output(command, shell=True)
		output = subprocess.check_output(command, shell=True)

		for i in [x.split(" ") for x in fq_stats.split("\n")[:-1]]:
			# Save fastq statistics
			# i[0] -> Filename
			save_eav(strain, "Number of Reads", i[1], Entity_Group = "Fastq Statistics", Sub_Entity = i[0], Tool="Awk")
			save_eav(strain, "Unique Reads", i[2], Entity_Group = "Fastq Statistics", Sub_Entity = i[0], Tool="Awk")
			save_eav(strain, "Frequency of Unique Reads", i[3], Entity_Group = "Fastq Statistics", Sub_Entity = i[0],  Tool="Awk")
			save_eav(strain, "Most abundant Sequence", i[4], Entity_Group = "Fastq Statistics",Sub_Entity = i[0],  Tool="Awk")
			save_eav(strain, "Number of times most abundant sequence occurs", i[5], Entity_Group = "Fastq Statistics", Sub_Entity = i[0],  Tool="Awk")
			save_eav(strain, "Frequency of Most Abundant Sequence", i[6], Entity_Group = "Fastq Statistics", Sub_Entity = i[0], Tool="Awk")

	for fq in fq_set:
		fastq_info = extract_fastq_info(fq)
		save_eav(strain, "Flowcell Lane",  fastq_info["flowcell_lane"], Entity_Group="Fastq_Meta", Sub_Entity= fq,  Tool="Python Script")
		save_eav(strain, "Index", fastq_info["index"], Entity_Group =  "Fastq_Meta",  Sub_Entity= fq,  Tool="Python Script")
		save_eav(strain, "Instrument", fastq_info["instrument"][1:], Entity_Group = "Fastq_Meta", Sub_Entity= fq, Tool="Python Script")
		save_eav(strain, "Read Pair", fastq_info["pair"], Entity_Group =  "Fastq_Meta",  Sub_Entity= fq, Tool="Python Script")

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
		reference_loc = "../reference/%s/%s.fa" % (reference, reference)
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
			save_eav(strain, i[1], i[2], Sub_Attribute=i[0], Entity_Group = "Bam Statistics", Sub_Entity = sample + ".bam", Tool = "Samtools Depth + Python")
		
		# Generate md5 of bam and store
		print bam_name
		save_md5([bam_name + ".bam"], Entity = strain)

		# Remove temporary files
		os.remove("%s.tmp.sorted.bam" % bam_name)
		os.remove("%s.tmp.bam" % bam_name)
		os.remove("%s.duplicate_report.txt" % bam_name)

		bam_set.append(bam_name + ".bam")

	# Merge BAMs into a single BAM.
	if len(bam_set) > 1:
		os.system("samtools merge -f -@ 4 %s.bam %s " % (strain, ' '.join(bam_set)))
		# Create record indicating the bams that were merged.
		save_eav(strain, sample + ".bam", ', '.join(bam_set), Sub_Entity = "Constitutive Bams", Entity_Group = "Bam Merged Statistics")
		save_eav(strain, sample + ".bam", len(bam_set), Sub_Entity = "Constitutive Bam Count", Entity_Group = "Bam Merged Statistics")
	else:
		os.system("mv %s %s.bam" % (bam_set[0] ,sample))
	os.system("samtools index %s.bam" % sample)

# Generate Depth and coverage statistics of the merged bam
if process_steps["bam_stats"] == True:
	for i in coverage(sample + ".bam", "chrM"):
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
		os.remove(b)
		os.remove(b + ".bai")
	except:
		pass


# Finally - Generate and Save an md5 sum.
save_md5([sample + ".bam"], Entity = strain)

# Move bam to bam folder.
os.system("mv %s.bam ../bam/%s.bam" % (strain, sample))
os.system("mv %s.bam.bai ../bam/%s.bam.bai" % (strain, sample))
os.chdir("../bam/") # Change Dirs into bam directory.

#====================#
# Variant Calling    #
#====================#

# Generate and store bcf stats
def gen_bcf_stats(strain, bcf, description, Sub_Attribute):
	# Generate stats - Following Het Polarization, Overwrites previous stat file.
	tmp =  tempfile.NamedTemporaryFile().name
	os.system("bcftools stats --samples - %s > %s" % (bcf, tmp))

	# Get variant statistics.
	for i in import_stats("%s" % (tmp)):
		save_eav(strain, i[0], i[1], Entity_Group = "BCF Statistics", Sub_Entity = "%s.%s.bcf" % (strain, description), Sub_Attribute = Sub_Attribute, Tool="bcftools")


# Fetch contigs
contigs = {}
for x in file("../reference/%s/%s.fa.fai" % (reference, reference), 'r').read().split("\n")[:-1]:
	contigs[x.split("\t")[0]] = int(x.split("\t")[1])

if process_steps["call_variants"] == True:
	command = "parallel --gnu 'samtools mpileup -t DP,DV,DP4,SP -g -f ../reference/%s/%s.fa -r {} %s.bam | bcftools call --format-fields GQ,GP -c -v > raw.%s.{}.bcf' ::: %s" % (reference, reference, strain, strain, ' '.join(contigs.keys()))
	print command
	os.system(command)
	os.system("echo %s | bcftools concat -O b `ls -v raw.%s.*.bcf` > %s.single.bcf" % (contigs.keys(), sample , sample.replace(".txt","")))
	
	# Index
	os.system("bcftools index -f %s.single.bcf" % sample)

	# Remove temporary files
	os.system("rm -f raw.%s.*" % sample)


gen_bcf_stats(strain, "%s.single.bcf" %  strain, "Single", "No Filter")

os.system("mv %s.single.bcf ../individual_bcf/%s.single.bcf" % (strain, strain))
os.system("bcftools index ../individual_bcf/%s.single.bcf" % (strain))
