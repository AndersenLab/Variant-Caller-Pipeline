#!/bin/bash
cd ../../data/ancillary/
wget --timestamping 'http://hgdownload.soe.ucsc.edu/goldenPath/ce10/database/rmsk.txt.gz' -O LCR_rmsk.txt.gz
gunzip -kfc LCR_rmsk.txt.gz