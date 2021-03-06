#!/bin/bash -v

# must invoke this script using bash 02b_radseq.sh

cd ../../data/radseq/

wget 'https://github.com/AndersenLab/Rad-Seq/blob/master/data/41188SNPset.txt?raw=true' -O 41188SNPset.txt

# Download files necessary to perform liftover.
wget ftp://ftp.sanger.ac.uk/pub2/wormbase/software/Remap-between-versions/remap.tar.bz2
gunzip remap.tar.bz2
tar -xf remap.tar

for r in `ls Remap-for-other-groups/*`; do
	echo $r
done;

mv Remap-for-other-groups/* .
rm -d Remap-for-other-groups/

# Convert to gff
awk 'NR>1{ split($1,a,"_"); print "CHROMOSOME_" substr($1,0,1) "\t.\t" $1 "\t" a[2] "\t" a[2] "\t.\t+\t.\t."}' 41188SNPset.txt | egrep -v "MtDNA" > 41188SNPset.ws210.gff
# liftover
perl remap_gff_between_releases.pl --release1=210 --release2=220 --gff=41188SNPset.ws210.gff --output=41188SNPset.ws220.gff

ws210=`mktemp -t tmp`
cut -f 3,4,5 41188SNPset.ws210.gff > $ws210

header=`mktemp -t tmp`
head -n 1 41188SNPset.txt | awk '{gsub(/^/,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"); print}' > $header

join -t $'\t' $ws210 41188SNPset.txt | gawk -f ../../scripts/prep/02b_radseq_awk.awk | cat ../ancillary/vcf_header_radseq.txt $header - | bcftools view -O b > andersen08_radseq.ws220.bcf


mv andersen08_radseq.ws220.bcf ../vcf/andersen08_radseq.ws220.bcf
bcftools index ../vcf/andersen08_radseq.ws220.bcf