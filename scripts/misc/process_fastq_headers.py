
from itertools import groupby as g
def most_common(L):
  return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]

f = filter(len,open("sum.txt").read().split("|"))[:-1]

q = open("BGI2_fastq_report.txt", "w")

for lines in filter(len,open("sum.txt").read().split("|"))[:-1]:
	line_set = filter(len,lines.split("\n"))
	file = line_set[0]
	lines = line_set[1:]

	# Pull out info
	machine = lines[1].split(':')[0].replace("@","")
	lane = lines[1].split(':')[2]
	index = most_common([x[x.index("#")+1:-2] for x in lines])
	
	q.write('\t'.join([file, machine, lane, index, "\n"]))
	