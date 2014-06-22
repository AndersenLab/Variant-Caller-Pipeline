#!/bin/bash

cd ../../data/

# Convert to gff
awk 'NR>1{ split($1,a,"_"); print "CHROMOSOME_" substr($1,0,1) "\t.\t" $1 "\t" a[2] "\t" a[2] "\t.\t+\t.\t."}' 41188SNPset.txt | egrep -v "MtDNA" > 41188SNPset.ws210.gff
# liftover
perl remap_gff_between_releases.pl --release1=210 --release2=225 --gff=41188SNPset.ws210.gff --output=41188SNPset.ws225.gff


join -t $'\t' <(cut -f 3,4,5 41188SNPset.ws220.gff) 41188SNPset.txt | \
gawk -f ../../scripts/prep/radseq_awk.awk | cat vcf_header.txt <(head -n 1 41188SNPset.txt | awk '{gsub(/^/,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"); print}') - | bgzip > andersen08_radseq.ws225.vcf.gz

mv andersen08_radseq.ws225.vcf.gz ../vcf/andersen08_radseq.ws225.vcf.gz
tabix ../vcf/andersen08_radseq.ws225.vcf.gz