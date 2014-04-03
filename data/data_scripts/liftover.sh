#!/bin/bash
# Download c. elegans genome assemblies ws220 and ws235
# Much thanks to http://genomewiki.ucsc.edu/images/9/91/SameSpeciesBlatSetup.sh.txt

export working_dir=../raw/assemblies

cd ${working_dir}

# Download ws220 (2bit file)
wget --timestamping 'http://hgdownload-test.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.2bit'
twoBitToFa ce10.2bit ce10.fa

# Download ws235 
wget --timestamping 'ftp://ftp.wormbase.org/pub/wormbase/releases/WS235/species/c_elegans/c_elegans.WS235.genomic.fa.gz'
gunzip c_elegans.WS235.genomic.fa.gz
mv c_elegans.WS235.genomic.fa ws235.fa
FaToTwoBit ws235.fa ws235.2bit

# Generate chromosome sizes
twoBitInfo ws235.2bit ws235.chrom.sizes
twoBitInfo ce10.2bit ce10.chrom.sizes

# Split the (new) genome
mkdir ws235 ce10
faSplit byname ws235.fa ws235/
faSplit byname ce10.fa ce10/

# Generate .ooc file; Contains a list of common -mers.
blat ce10.2bit /dev/null /dev/null -tileSize=11 -makeOoc=11.ooc -repMatch=100


############ 4/3/14 END #############

blat ${working_dir}/ce10.2bit ${working_dir}/ws220.2bit OLD.chr0.psl -tileSize=11 -minScore=100 -minIdentity=98 

#######################################################################
# First part of the procedure is a blat run between target and query chunks
# decide on a work-directory where this will take place:
export workDirectory="../raw/temp"
# two directories will be created here where work-steps take place:

mkdir -p ${workDirectory}/run.blat ${workDirectory}/run.chain

# the resulting liftOver chain file will end up in ${workDirectory}/

# location where your two .2bit files are located:
export twoBitDirectory="../../assemblies"

#######################################################################
# targetDb is the sequence you have good coordinates you want to convert
export targetDb="ce10"
# typical chunk size for a target sequence is 10,000,000 bases
export targetChunkSize="10000000"
# no need to overlap chunks, either target or query:
export chunkOverlapSize="0"
# allow only a single sequence to exist in a single chunk group
export chunkCountLimit="1"

# queryDb is the sequence you want to convert coordinates from the target
#    to the query
export queryDb="ws220"
# typical chunk size for a query sequence is also 10,000,000 bases
export queryChunkSize="10000000"

#######################################################################
# working in the run.blat directory; Generate chrom.sizes

cd ${workDirectory}/run.blat

# directory tParts will be created, make sure it is clean
rm -fr tParts
twoBitInfo ${twoBitDirectory}/${targetDb}.2bit stdout | sort -k2nr > ${targetDb}.chrom.sizes
twoBitInfo ${twoBitDirectory}/${queryDb}.2bit stdout | sort -k2nr > ${queryDb}.chrom.sizes

#######################################################################
# Construct the 11.ooc file for the target sequence, assuming

#  the printout of the faSize is a line something like:
# 23332831 bases (10 N's 23332821 real 18400480 upper 4932341 lower) in 16 sequences in 1 files
# use the fifth field there, the 23332821 "real" number of bases:
export realSize=`twoBitToFa ${twoBitDirectory}/${targetDb}.2bit stdout | faSize stdin | grep real | awk '{print $5}'`
# Use a repMatch of 100:


blat $(path1) $(path2) {check out line+ ../raw/$(root1)_$(root2).psl} -tileSize=12 -ooc=11.ooc -minScore=100 -minIdentity=98 -fastMap

#######################################################################
# Construct the result output PSL directories

cd ${workDirectory}/run.blat
mkdir ${workDirectory}/run.blat/psl
for F in `cat t.lst`
do
  B=`basename ${F}
  mkdir ${workDirectory}/run.blat/psl/${B}
done

#######################################################################
# construct gensub2 template file
#  ( allow only ${workDirectory} to shell expand here, the other $ variables
#	are for gensub2 use, hence \$ to keep them literal )
# the "blatJob.csh" script can be found at:

cat << EOF > template
#LOOP
./blatJob.csh \$(path1) \$(path2) ${workDirectory}/run.blat/psl/\$(file1)/\$(file2).psl
#ENDLOOP
EOF

## construct jobList for each query chunk blat to each target chunk:

gensub2 t.lst q.lst template jobList

# the resulting jobList is a listing of commands to be run which will perform
# the blat on each specified target/query chunk pair of sequences
# With the blatJob.csh in place in this run.blat directory, you can simply
# run each job in the jobList:
#    chmod +x ./jobList
#    ./jobList > jobs.log 2>&1

#######################################################################