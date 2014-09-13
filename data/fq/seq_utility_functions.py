import gzip, re, subprocess
from peewee import *
from itertools import groupby as g


#===========#
# Functions #
#===========#


def most_common(L):
  return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]

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

#print import_stats("JU1395.tmp.stats.txt")

def save_md5(files = [], type = ""):
	md5 = subprocess.check_output("parallel %s ::: %s" % (md5_system[system_type], ' '.join(files)), shell=True)
	m = re.findall('MD5 \((.*)\) = (.*)', md5)
	for i in m:
		# Check File Hashes
		store_eav("MD5 Hash", i[0], "File Hash", i[1], Tool="MD5")