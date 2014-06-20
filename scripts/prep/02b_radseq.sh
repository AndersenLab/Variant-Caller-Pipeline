#!/bin/bash

cd ../../data/

# Convert to gff
awk 'NR>1{ split($1,a,"_"); print "CHROMOSOME_" substr($1,0,1) "\t.\t" $1 "\t" a[2] "\t" a[2] "\t.\t+\t.\t."}' 41188SNPset.txt | egrep -v "MtDNA" > 41188SNPset.ws210.gff
# liftover
perl remap_gff_between_releases.pl --release1=210 --release2=220 --gff=41188SNPset.ws210.gff --output=41188SNPset.ws220.gff


join -t $'\t' <(cut -f 3,4,5 41188SNPset.ws220.gff) 41188SNPset.txt | awk '{ ORS=""; for(item in freq) freq[item]=0; for ( i=4; i<NF; i++ ) freq[$i]++; for(item in freq) print "\titem\t" item "\tfreq\t" freq[item]; print "\n"; print "chr" gensub("_.*$","","g",$1) "\t" $2 "\torig_" $1 "\txREF\txALT\t.\t.\t.\t.\t"; for ( i=4; i<NF; i++ ) print $i "\t"; print "\n"}' |cat <(head -n 1 41188SNPset.txt | awk '{gsub(/^/,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"); print}') - | head -n 20

awk '{print "chr" gensub("_.*$","","g",$3) "\t" $4 "\t" $3}' 41188SNPset.ws220.gff 