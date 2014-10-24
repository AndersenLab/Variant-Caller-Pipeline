#!/bin/bash

# This script can be used to compare andersen sequencing with mmp.

# Download mmp data.
cd ../../data/mmp/
cp /Volumes/PortusTutus/backup/mmp_isolate_vcf.tar.gz isolate_vcf.tar.gz
gunzip isolate_vcf.tar.gz
tar -xf isolate_vcf.tar 
rm isolate_vcf.tar

# Fix ED3017 - appears to be a duplicate line or something weird.
sed '107577d' ED3017.wb225.fa.combined.vcf | sed '107576d' ED3017.wb225.fa.combined.vcf > ED3017.f.vcf
rm ED3017.wb225.fa.combined.vcf
mv ED3017.f.vcf ED3017.wb225.fa.combined.vcf


function fix_vcf (){
    STRAIN_NAME=${1%%\.*}
    awk -v strain=$STRAIN_NAME -F $'\t' 'BEGIN {OFS = FS}
     $1 ~ "^\##" { print }
     $1 ~ "^\#CHROM" { print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; print $0 "\tFORMAT\t" strain}
     $1 !~ "^\#" { $3="."; $8 = "."; print "chr" $0 "\tGT\t1/1"}
     ' $1 | sed -e 's/tDNA//g' > $STRAIN_NAME.fixed.vcf
    bgzip -f $STRAIN_NAME.fixed.vcf
    tabix -f $STRAIN_NAME.fixed.vcf.gz
}
export -f fix_vcf

parallel fix_vcf ::: `ls *.vcf`

# Combine vcfs; collapse variants.
vcf-merge -R 0/0 -c both `ls *.fixed.vcf.gz`  > mmp.tmp.vcf


# Generate a new, fixed header

egrep '^#' mmp.tmp.vcf | sed -e 's/contig=<ID=/contig=<ID=chr/g' | sed -e 's/tDNA//g' | \
sed 's/ED3021/ED3005/g' | \
sed 's/ED3042/ED3048/g' | \
sed 's/JU1401/JD1409/g' | \
sed 's/JU263/JU310/g' | \
sed 's/JU322/JU323/g' | \
sed 's/JU361/JU367/g' | \
sed 's/JU533/JU1213/g' | \
sed 's/MY14/MY23/g' | \
sed 's/MY6/MY18/g' | \
sed 's/PX174/RC301/g' > corrected_header.txt


# Fix Header and combine
cat corrected_header.txt <(grep '#' -v mmp.tmp.vcf) > mmp.vcf

# Compres and Index.
bgzip mmp.vcf
tabix mmp.vcf.gz

mv mmp.vcf.gz ../vcf/mmp.vcf.gz
mv mmp.vcf.gz.tbi ../vcf/mmp.vcf.gz.tbi

