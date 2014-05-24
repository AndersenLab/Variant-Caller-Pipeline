#!/usr/bin/python
import re
from itertools import groupby as g
import subprocess

def most_common(L):
  return max(g(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]

r = subprocess.check_output("""
for r in `ls *.gz`; 
do
echo "$r" 
gunzip -cq $r | head -n 1000 | grep '^@' - 
echo "|" 
done;
""", shell=True)

prog = re.compile('@.*:')

f = open("fastq_summary.txt", "w")

f.write("File\tMachine\tLane\tIndex\tPair\n")

for i in [filter(len,x.split('\n')) for x in r.split("|")][:-1]:
	f.write('\t'.join([i[0], re.match("@(.*?):", i[1]).group(1), re.match("@(.*?):([0-9]+)", i[1]).group(2), most_common([x[-10:-2] for x in i[1:]]),  i[1][-1:], "\n"]))