#!/bin/bash

# Download mmp data.
cd ../../data/mmp/
wget --timestamping 'http://waterston.gs.washington.edu/trackhubs/data/isolate_vcf.tar.gz'
gunzip isolate_vcf.tar.gz
tar -xf isolate_vcf.tar 
rm isolate_vcf.tar

# Add ALT genotype column 1/1 to end; 
for r in `ls *.vcf`; do
    echo $r
    STRAIN_NAME=${r%%\.*}
    awk -v strain=$STRAIN_NAME '
     $1 ~ "^\##" { print }
     $1 ~ "^\#CHROM" { print $0 "\t" strain}
     $1 !~ "^\#" { $3="."; $8="."; print $0 "\t1|1"}
     ' $r > $STRAIN_NAME.fixed.vcf
    bgzip -f $STRAIN_NAME.fixed.vcf
    tabix $STRAIN_NAME.fixed.vcf.gz
done
# Combine vcfs; collapse variants.
vcf-merge -c both `ls *.fixed.vcf.gz` > mmp.vcf

# Compres and Index.
bgzip mmp.vcf
tabix mmp.vcf.gz

mv mmp.vcf.gz ../vcf/mmp.vcf.gz
mv mmp.vcf.gz.tbi ../vcf/mmp.vcf.gz.tbi

