!#/bin/bash
# This script can be used to extract information from fastqs.

for r in `ls *.gz`; 
do 
echo "$r"  >> sum.txt
gunzip -c $r | head -n 100 | egrep '^@'   >> sum.txt; 
echo "|" >> sum.txt
done;
