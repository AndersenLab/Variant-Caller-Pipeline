#!/bin/bash

cd /Volumes/PortusTutus/Raw/SRA_FASTQ_MMP/fq

# Fetch Run Numbers

function SRX_fetch_fastq() {
    sra_set=`esearch -db sra -query $1 | efetch -format docsum | xtract -element Run@acc`
    echo "Downloading Run $1:" 
    echo ${sra_set}
    echo "-------"
    for SRA in $sra_set; do
        echo "Downloading $SRA"
        fastq-dump $SRA
    done;
}


for r in `cat MMP_ids.txt`; do
    SRX_fetch_fastq $r
done;