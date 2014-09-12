import gzip, re
from itertools import groupby as g

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