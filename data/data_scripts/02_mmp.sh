#!/bin/bash

cd ../raw/mmp

# Download Data from "The million mutation project: A new approach to genetics in Caenorhabditis elegans"
wget 'http://genome.sfu.ca/mmp/mmp_wild_isolate_data_Mar13.txt'

head -n 200 mmp235.gff > mmp_test.gff

# Convert file to a 'gff'
awk -F $'\t' '{ print $4 "\t.\t.\t" $5 "\t" ($5+1) "\t.\t-\t" $1 }'  <(tail +2 mmp_wild_isolate_data_Mar13.txt) > mmp235.gff


# File needs to be in gff format.

# Download the necessary files for the conversion
# Thanks to https://github.com/pruzanov/ME.scripts/blob/master/remap_gff/unmap_gff_between_releases.pl

perl remap_gff_between_releases.pl -gff mmp235.gff -output mmp220.gff  -release1 195 -release2 235

# Check if they are the same...
md5 mmp220.gff
md5 mmp235.gff