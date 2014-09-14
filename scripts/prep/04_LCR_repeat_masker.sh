#!/bin/bash
## LCR 1
cd ../../data/LCR_Region/
wget 'http://hgdownload.soe.ucsc.edu/goldenPath/ce10/database/rmsk.txt.gz' -O LCR_rmsk.txt.gz
gunzip -kfc LCR_rmsk.txt.gz | grep 'Low_complexity' | cut -f 6,7,8 > LCR_temp.bed
rm LCR_rmsk.txt.gz

# Download c. elegans chromosome information
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
	"select chrom, size from ce10.chromInfo"  > ce10.genome

bedtools complement -i LCR_temp.bed  -g ce10.genome | sort -k 1,1 -k2,2n > ce10.ucsc.masked.bed

## LCR 2

wget 'ftp://ftp.wormbase.org/pub/wormbase/releases/WS220/species/c_elegans/c_elegans.WS220.genomic_masked.fa.gz' 
python ../../scripts/misc/generate_masked_ranges.py c_elegans.WS220.genomic_masked.fa.gz CHROMOSOME_ chr | sed 's/chrMtDNA/chrM/g' > WS220.tmp.bed

bedtools complement -i WS220.tmp.bed  -g ce10.genome | sort -k 1,1 -k2,2n > WS220.wormbase.masked.bed

## Generate Stats on repeat masking, and store in database.
bedtools genomecov -i ce10.ucsc.masked.bed -g ce10.genome > ce10.ucsc.masked.genomecov.statsitics.txt
bedtools genomecov -i WS220.wormbase.masked.bed -g ce10.genome > WS220.wormbase.masked.genomecov.statistics.txt



rm c_elegans.WS220.genomic_masked.fa.gz
rm LCR_temp.bed
rm WS220.tmp.bed
