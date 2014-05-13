#!/bin/bash

# Download mmp data.
cd ../../data/mmp/
wget --timestamping 'http://waterston.gs.washington.edu/trackhubs/data/isolate_vcf.tar.gz'
gunzip isolate_vcf.tar.gz
tar -xf isolate_vcf.tar 
rm isolate_vcf.tar

# Fix ED3017 - appears to be a duplicate line or something.
sed '107577d' ED3017.wb225.fa.combined.vcf | sed '107576d' ED3017.wb225.fa.combined.vcf > ED3017.f.vcf
rm ED3017.wb225.fa.combined.vcf
mv ED3017.f.vcf ED3017.wb225.fa.combined.vcf

# Add ALT genotype column 1/1 to end; 
for r in `ls *.vcf`; do
    echo $r
    STRAIN_NAME=${r%%\.*}
    awk -v strain=$STRAIN_NAME -F $'\t' 'BEGIN {OFS = FS}
     $1 ~ "^\##" { print }
     $1 ~ "^\#CHROM" { print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; print $0 "\tFORMAT\t" strain}
     $1 !~ "^\#" { $3="1"; $8 = "."; print $0 "\tGT\t1|1"}
     ' $r > $STRAIN_NAME.fixed.vcf
    bgzip -f $STRAIN_NAME.fixed.vcf
    tabix $STRAIN_NAME.fixed.vcf.gz
done
# Combine vcfs; collapse variants.
vcf-merge -R 0/0 -c both `ls *.fixed.vcf.gz` > mmp.vcf

# Compres and Index.
bgzip mmp.vcf
tabix mmp.vcf.gz

mv mmp.vcf.gz ../vcf/mmp.vcf.gz
mv mmp.vcf.gz.tbi ../vcf/mmp.vcf.gz.tbi

