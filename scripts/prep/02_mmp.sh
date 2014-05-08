#!/bin/bash

# Download mmp data.
cd ../../data/mmp/
wget --timestamping 'http://waterston.gs.washington.edu/trackhubs/data/isolate_vcf.tar.gz'
gunzip isolate_vcf.tar.gz
tar -xf isolate_vcf.tar 
rm isolate_vcf.tar

# bgzip
for r in `ls *.vcf`; do
	bgzip $r
	tabix $r.gz
done
# Combine vcfs
vcf-merge `ls *.vcf.gz` > mmp.vcf

# Remove extra files. Only keep mmp.vcf
rm `find . -type f -not -name "mmp.vcf"`