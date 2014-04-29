!#/bin/bash
# Use to extract info from a set of fastq's

for r in `ls *.gz`; 
do 
echo "$r"  >> sum.txt
gunzip -c $r | head -n 100 | egrep '^@'   >> sum.txt; 
echo "|" >> sum.txt
done;
