#!/bin/bash -v

# must invoke this script using bash 02b_radseq.sh

cd ../../data/radseq/

# Download files necessary to perform liftover.
wget ftp://ftp.sanger.ac.uk/pub2/wormbase/software/Remap-between-versions/remap.tar.bz2
gunzip remap.tar.bz2
tar -xf remap.tar

for r in `ls Remap-for-other-groups/*`; do
	echo $r
done;

mv Remap-for-other-groups/* .
rm -d Remap-between-versions/

# Convert to gff
awk 'NR>1{ split($1,a,"_"); print "CHROMOSOME_" substr($1,0,1) "\t.\t" $1 "\t" a[2] "\t" a[2] "\t.\t+\t.\t."}' 41188SNPset.txt | egrep -v "MtDNA" > 41188SNPset.ws210.gff
# liftover
perl remap_gff_between_releases.pl --release1=210 --release2=225 --gff=41188SNPset.ws210.gff --output=41188SNPset.ws225.gff

ws210=`mktemp -t tmp`
cut -f 3,4,5 41188SNPset.ws210.gff > $ws210

header=`mktemp -t tmp`
head -n 1 41188SNPset.txt | awk '{gsub(/^/,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"); print}' > $header

join -t $'\t' $ws210 41188SNPset.txt | gawk -f ../../scripts/prep/02b_radseq_awk.awk | cat vcf_header.txt $header - | bgzip > andersen08_radseq.ws225.vcf.gz

mv andersen08_radseq.ws225.vcf.gz ../vcf/andersen08_radseq.ws225.vcf.gz
tabix ../vcf/andersen08_radseq.ws225.vcf.gz