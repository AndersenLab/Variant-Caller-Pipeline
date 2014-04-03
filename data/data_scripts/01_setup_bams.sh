#!/bin/bash
# Download BGI1, BGI2, BGI3 + Princeton

# Setup Dirs
bamdir="../raw/bam"
mkdir "$bamdir/bam_headers"

wget --timestamping --directory-prefix $bamdir -i ../raw/bam_urls/BGI1_urls.txt
wget --timestamping --directory-prefix $bamdir -i ../raw/bam_urls/BGI2_urls.txt
wget --timestamping --directory-prefix $bamdir -i ../raw/bam_urls/BGI3_Princeton_urls.txt

# Rename appropriately; Output BAM headers
for file in `ls $bamdir/*.bam*`
do
    if [ -f $FILE ]; then
      # create a variable that strips URL params from end of filename.
      sfile=${file/?ticket=t_euEWjVvI/}
      mv "$file" "$sfile"
      # Output sam headers and rename extension: 
      samtools view -H $sfile > "$bamdir/bam_headers/$(basename ${sfile/.bam/.bam_header})"
    fi
done

# Check that bam filenames against bam header information.
for f in `ls $bamdir/bam_headers/*`
    do
        # Compare header to filename
        header_spec=$(egrep "SM:[A-Za-z_0-9\.]+" --only-matching  $f)
        if [[ "${header_spec/SM:/}" != "$(basename ${f/.bam_header/})" ]]
         then echo "Mismatch header:$header_spec filename:$(basename ${f/.bam_header/})"
         else echo "Okay"
        fi
done

# Download and install liftover
wget -O "liftover" http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/liftOver
mv liftover /usr/local/bin

# Download Data from "The million mutation project: A new approach to genetics in Caenorhabditis elegans"
wget 