#!/bin/bash
# Download c. elegans genome assemblies ws220 and ws235
# Much thanks to http://genomewiki.ucsc.edu/images/9/91/SameSpeciesBlatSetup.sh.txt


#---------------------------------------------------------------------------#
# First part of the procedure is a blat run between target and query chunks #
# decide on a work-directory where this will take place:                    #
#---------------------------------------------------------------------------#
export workDirectory="raw/liftover"
mkdir ${workDirectory}
# two directories will be created here where work-steps take place:
mkdir ${workDirectory}/run.blat ${workDirectory}/run.chain



#----------#
# Download #
#----------#
cd ${workDirectory}
# Download ws220 (2bit file)
wget --timestamping --directory-prefix "assemblies" 'http://hgdownload-test.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.2bit'

# Download ws235 
wget --timestamping --directory-prefix "assemblies"  'ftp://ftp.wormbase.org/pub/wormbase/releases/WS235/species/c_elegans/c_elegans.WS235.genomic.fa.gz'
gunzip assemblies/c_elegans.WS235.genomic.fa.gz
mv assemblies/c_elegans.WS235.genomic.fa assemblies/ws235.fa
FaToTwoBit assemblies/ws235.fa assemblies/ws235.2bit


#------------------------------------------------------------------------#
# targetDb is the sequence you have good coordinates you want to convert #
#------------------------------------------------------------------------#
export targetDb="ws235"
# typical chunk size for a target sequence is 10,000,000 bases
export targetChunkSize="10000000"
# no need to overlap chunks, either target or query:
export chunkOverlapSize="0"
# allow only a single sequence to exist in a single chunk group
export chunkCountLimit="1"

# queryDb is the sequence you want to convert coordinates from the target
#    to the query
export queryDb="ce10"
# typical chunk size for a query sequence is also 10,000,000 bases
export queryChunkSize="10000000"

#------------------------------------------------------------------------#
# Compute partitions (coordinate ranges) for cluster jobs.               #
# working in the run.blat directory                                      #
#------------------------------------------------------------------------#

cd run.blat

# location where your two .2bit files are located:
export twoBitDirectory="../assemblies"

# directory tParts will be created, make sure it is clean
rm -fr tParts
twoBitInfo ${twoBitDirectory}/${targetDb}.2bit stdout | sort -k2nr > ${targetDb}.chrom.sizes
# check the usage message for this command in the kent source tree
#   to see what all the arguments mean:
../../../data_scripts/partitionSequence.pl ${targetChunkSize} \
   ${chunkOverlapSize} ${twoBitDirectory}/${targetDb}.2bit ${targetDb}.chrom.sizes \
   ${chunkCountLimit} \
   -lstDir=tParts > t.lst
# output t.lst is the listing of the target chunks to work with

# same procedure for the query bits:
rm -fr qParts
twoBitInfo ${twoBitDirectory}/${queryDb}.2bit stdout | sort -k2nr > ${queryDb}.chrom.sizes
# check the usage message for this command in the kent source tree
#   to see what all the arguments mean:
../../../data_scripts/partitionSequence.pl ${queryChunkSize} \
   ${chunkOverlapSize} ${twoBitDirectory}/${queryDb}.2bit ${queryDb}.chrom.sizes \
   ${chunkCountLimit} \
   -lstDir=qParts > q.lst
# output q.lst is the listing of the query chunks to work with

#######################################################################
# Construct the 11.ooc file for the target sequence, assuming

#  the printout of the faSize is a line something like:
# 23332831 bases (10 N's 23332821 real 18400480 upper 4932341 lower) in 16 sequences in 1 files
# use the fifth field there, the 23332821 "real" number of bases:
export realSize=`twoBitToFa ${twoBitDirectory}/${targetDb}.2bit stdout | faSize stdin | grep real | awk '{print $5}'`
# calculate repMatch based on a proportion of the hg19 genoe size
#export repMatch=`calc \( ${realSize} / 2897316137 \) \* 1024`
# That gives a really small number for small genomes, ignore that result
# and simply use a repMatch of 100:
blat ${twoBitDirectory}/${targetDb}.2bit /dev/null \
        /dev/null -tileSize=11 -makeOoc=11.ooc -repMatch=100

#######################################################################
# Construct the result output PSL directories

mkdir psl
for F in `cat t.lst`;
do
  B=`basename ${F}`
  mkdir psl/${B}
done

#######################################################################
# construct gensub2 template file
#  ( allow only ${workDirectory} to shell expand here, the other $ variables
#	are for gensub2 use, hence \$ to keep them literal )
# the "blatJob.csh" script can be found at:

export twoBitDirectory="../assemblies"

cat << EOF > template
#LOOP
./blatJob.csh \$(path1) \$(path2) ${workDirectory}/run.blat/psl/\$(file1)/\$(file2).psl
#ENDLOOP
EOF

## construct jobList for each query chunk blat to each target chunk:



gensub2 t.lst q.lst template jobList
chmod +x jobList 
# the resulting jobList is a listing of commands to be run which will perform
# the blat on each specified target/query chunk pair of sequences
# With the blatJob.csh in place in this run.blat directory, you can simply
# run each job in the jobList:
#    chmod +x ./jobList
#    ./jobList > jobs.log 2>&1

#######################################################################