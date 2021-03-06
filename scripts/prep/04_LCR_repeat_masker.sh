#!/bin/bash
## LCR 1
#cd ../../data/LCR_Region/
#curl 'http://hgdownload.soe.ucsc.edu/goldenPath/ce10/database/rmsk.txt.gz' | gunzip -kfc - | grep 'Low_complexity' | cut -f 6,7,8 | bedtools merge -i - | awk '{ print $0 "\tLCR_UCSC_ce10_ws220"}'  > ce10.ucsc.masked.bed

# Download c. elegans chromosome information
#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
#	"select chrom, size from ce10.chromInfo"  > ce10.genome

## LCR 2

curl 'ftp://ftp.wormbase.org/pub/wormbase/releases/WS235/species/c_elegans/c_elegans.WS235.genomic_masked.fa.gz' > c_elegans.WS235.genomic_masked.fa.gz
python ../../scripts/misc/generate_masked_ranges.py c_elegans.WS235.genomic_masked.fa.gz | awk '{ print $0 "\tLCR_Wormbase_235"}' >  WS235.wormbase.masked.bed

## Generate Stats on repeat masking, and store in database.
#bedtools genomecov -i ce10.ucsc.masked.bed -g ce10.genome > ce10.ucsc.masked.genomecov.statsitics.txt
bedtools genomecov -i WS235.wormbase.masked.bed -g ws235.genome > WS235.wormbase.masked.genomecov.statistics.txt

#bgzip -f ce10.ucsc.masked.bed
bgzip -f WS235.wormbase.masked.bed
#tabix -f ce10.ucsc.masked.bed.gz
tabix -f WS235.wormbase.masked.bed.gz

#rm c_elegans.WS220.genomic_masked.fa.gz
#rm LCR_temp.bed
#rm WS220.tmp.bed
