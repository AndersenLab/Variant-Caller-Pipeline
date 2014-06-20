#!/bin/bash
cd ../../data/ancillary/
wget 'http://hgdownload.soe.ucsc.edu/goldenPath/ce10/database/rmsk.txt.gz' -O LCR_rmsk.txt.gz
gunzip -kfc LCR_rmsk.txt.gz | cut -f 6,7,8 > LCR_ce10_rmsk.bed
rm LCR_rmsk.txt.gz