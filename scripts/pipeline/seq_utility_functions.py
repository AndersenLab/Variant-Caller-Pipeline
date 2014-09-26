import gzip, re, subprocess, os
from peewee import *
from itertools import groupby as g


#===========#
# Functions #
#===========#


def most_common(L):
  try:
  	return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]
  except:
  	return ""


def extract_fastq_info(fastq):
	"""
	This function will extract information from the header lines of a demultiplexed fastq. Requires gzip to be imported.
	"""
	f = gzip.open(fastq, 'rb')
	header_lines = [x.replace("\n","") for x in f.readlines(10000) if x.startswith("@")]

	for heading in header_lines:
			l = re.split(r'(\:|#| )',heading)
			line = {}
			index_set = []
			if len(l) == 11:
				line["instrument"] = l[0]
				line["flowcell_lane"] = l[2]
				line["flowcell_tile"] = l[4]
				line["x_coord"] = l[6]
				line["y_coord"] = l[8]
				try:
					line["pair"] = l[10].split("/")[1]
					index_set.append(l[10].split("/")[0])
				except:
					pass
			elif len(l) == 21:
				line["instrument"] = l[0]
				line["run_id"] = l[2]
				line["flowcell_id"] = l[4]
				line["flowcell_lane"] = l[6]
				line["flowcell_tile"] = l[8]
				line["x_coord"] = l[10]
				line["y_coord"] = l[12]
				line["pair"] = l[14]
				line["filtered"] = l[16]
				line["control_bits"] = l[16]
				line["index"] = l[20]
				index_set.append(l[20])
			else:
				print "error", l
			line["index"] = most_common(index_set)
	return line


import re

def import_stats(stat_file):
	"""
		Single Sample BCF
		Simple script for importing some statistics from a bcf file - mostly discrete.
	"""
	f = file(stat_file, 'r')
	data = []
	single_stat_sets = ["PSC", "TSTV"]
	with open(stat_file, 'r+') as f:
		for line in f:
			line = line.replace("\n", "")
			if line.startswith("SN"):
				SN_Data = line.replace(":","").split("\t")[2:4]
				data += [(SN_Data[0].title(), SN_Data[1])]
			elif any([line.startswith("# %s\t" % x) for x in single_stat_sets]):
				key_group = line.split("\t")[0].replace("# ","")
				keys = ["%s - %s" % (key_group, x) for x in re.sub(r"\[[0-9]+\]", "", line).title().split("\t")[2:]]
				vals = f.next().replace("\n","").split("\t")[2:]
				data += zip(keys,vals)
	return data


import os
import re
import subprocess

def coverage(bam, mtchr = None):
	# Check to see if file exists
	if os.path.isfile(bam) == False:
		raise Exception("Bam file does not exist")
	header = subprocess.check_output("samtools view -H %s" % bam, shell = True)
	# Extract contigs from header and convert contigs to integers
	contigs = {}
	for x in re.findall("@SQ	SN:(?P<chrom>[A-Za-z0-9]+)\WLN:(?P<length>[0-9]+)", header):
		contigs[x[0]] = int(x[1])
		# Calculate Coverage for each chromosome individually
	coverage_dict = {}
	for c in contigs.keys():
		command = "samtools depth -r %s %s | awk '{sum+=$3;cnt++}END{print cnt \"\t\" sum}'" % (c, bam)
		coverage_dict[c] = {}
		print subprocess.check_output(command, shell = True).strip().split("\t")
		coverage_dict[c]["Bases Mapped"], coverage_dict[c]["Sum of Depths"] = map(int,subprocess.check_output(command, shell = True).strip().split("\t"))
		coverage_dict[c]["Breadth of Coverage"] = coverage_dict[c]["Bases Mapped"] / float(contigs[c])
		coverage_dict[c]["Depth of Coverage"] = coverage_dict[c]["Sum of Depths"] / float(contigs[c])
		coverage_dict[c]["Length"] = int(contigs[c])

	# Calculate Genome Wide Breadth of Coverage and Depth of Coverage
	genome_length = float(sum(contigs.values()))
	coverage_dict["genome"] = {}
	coverage_dict["genome"]["Length"] = int(genome_length)
	coverage_dict["genome"]["Bases Mapped"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k != "genome"])
	coverage_dict["genome"]["Sum of Depths"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k != "genome"])
	coverage_dict["genome"]["Breadth of Coverage"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k != "genome"]) / genome_length
	coverage_dict["genome"]["Depth of Coverage"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k != "genome"]) / genome_length

	if mtchr != None:
		# Calculate nuclear breadth of coverage and depth of coverage
		ignore_contigs = [mtchr, "genome", "nuclear"]
		coverage_dict["nuclear"] = {}
		coverage_dict["nuclear"]["Length"] = sum([x["Length"] for k,x in coverage_dict.iteritems() if k not in ignore_contigs ])
		coverage_dict["nuclear"]["Bases Mapped"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs])
		coverage_dict["nuclear"]["Sum of Depths"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs])
		coverage_dict["nuclear"]["Breadth of Coverage"] = sum([x["Bases Mapped"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs]) / genome_length
		coverage_dict["nuclear"]["Depth of Coverage"] = sum([x["Sum of Depths"] for k, x in coverage_dict.iteritems() if k not in ignore_contigs]) / genome_length

		# Calculate the ratio of mtDNA depth to nuclear depth
		mt_ratio = coverage_dict[mtchr]["Sum of Depths"] / float(coverage_dict["nuclear"]["Depth of Coverage"])

	# Flatten Dictionary 
	coverage = []
	for k,v in coverage_dict.items():
		for x in v.items():
			coverage += [(k,x[0], x[1])]
	return coverage



def save_md5(files = [], type = ""):
	md5 = subprocess.check_output("parallel %s ::: %s" % (md5_system[system_type], ' '.join(files)), shell=True)
	m = re.findall('MD5 \((.*)\) = (.*)', md5)
	for i in m:
		# Check File Hashes
		store_eav("MD5 Hash", i[0], "File Hash", i[1], Tool="MD5")