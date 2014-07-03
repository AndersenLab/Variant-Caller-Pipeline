#!/bin/bash
cd ../../data/ancillary/
wget 'http://hgdownload.soe.ucsc.edu/goldenPath/ce10/database/rmsk.txt.gz' -O LCR_rmsk.txt.gz
gunzip -kfc LCR_rmsk.txt.gz | grep 'Low_complexity' | cut -f 6,7,8 > LCR_ce10_rmsk.bed
rm LCR_rmsk.txt.gz

# Download c. elegans chromosome information
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
	"select chrom, size from ce10.chromInfo"  > ce10.genome

bedtools complement -i LCR_ce10_rmsk.bed  -g ce10.genome | sort -k 1,1 -k2,2n > ce10_no_LCR.bed