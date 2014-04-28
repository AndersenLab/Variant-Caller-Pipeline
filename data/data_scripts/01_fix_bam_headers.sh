#!/bin/bash
# Download BGI1, BGI2, BGI3 + Princeton

# Setup Dirs
bamdir="../raw/bam"
mkdir "$bamdir/bam_headers"

# Retrieve Bams, rename appropriately
sh 01a_download_rename_bams.sh


# Pull out headers
for file in `ls $bamdir/*.bam*`
do
    if [ -f $FILE ]; then
      # Output sam headers and rename extension: 
      samtools view -H $file > "$bamdir/bam_headers/$(basename ${file/.bam/.bam_header})"
    fi
done

# Check that bam filenames against bam header information.
for f in `ls $bamdir/bam_headers/*`
    do
        # Compare header to filename
        header_spec=$(egrep "SM:[A-Za-z_0-9\.]+" --only-matching  $f)
        if [[ "${header_spec/SM:/}" != "$(basename ${f/.bam_header/})" ]]
         then echo "Mismatch header:$header_spec filename:$(basename ${f/.bam_header/})"
        fi
done

# Download and install UCSC utilities.
for p in liftOver faSplit liftUp faToTwoBit twoBitToFa fetchChromSizes faSize
do
  wget -O "$p" "http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/$p"
  chmod 775 "$p"
  mv "$p" "/usr/local/bin"
done
