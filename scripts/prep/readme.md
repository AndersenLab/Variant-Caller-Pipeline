# Files in Prep/

* __00_demultiplex.sh__ - This script illustrates how seqtk was used to demultiplex the Princeton Data
* __00_prep_mmp.sh__ - Download million mutations project (MMP) data from SRA, convert to fastq, quality control, and rename in andersen lab Format.
* __Rockman_rename.sh__ - Script for renmaing Matt Rockman fastqs to Andersen lab format.
* __01_cleanup_fastqs.sh__ Renaming files from BGI, including fixing index errors.
* __02a_mmp_vcf.sh__ - Used to generate a vcf of Million Mutations Project (MMP) Data.
* __02b_radseq_awk.awk__ - Awk script used for processing radseq data.
* __02b_radseq.sh__ - Shell script for converting a CSV of Andersen lab Radseq data into a vcf for comparison with other vcfs.
* __03_construct_bam_sets.py__ - Script for constructing bam sets used in joint calling.
* __04_LCR_repeat_masker.sh__ - Script used to generate LCR files.
* __renaming_mmp__ - Excel file detailing how sequencing libraries were derived from fastq headers from SRA data in the Million Mutations Project.
