# Create a file with the princeton barcodes

cd "/Volumes/PortusTutus/Raw/FASTQ_clean/Princeton"

echo  """AF16	AGGTACCA
CB4856_CGC	ACCAACTA
DL238	AGTCAGAA
HK104	TGACGTCA
JU258	AATTCAAA
MY23	CTGCGACA
N2_CGC	TCGCAGGA
QX1430	CATCCGGA""" > princeton_barcodes.txt

a1="andersen_library1-for-210-cycles-d1df6acxx_6_read_1_passed_filter.fastq.gz"
a2="andersen_library1-for-210-cycles-d1df6acxx_6_read_3_passed_filter.fastq.gz"
aindex="andersen_library1-for-210-cycles-d1df6acxx_6_read_2_index_read_passed_filter.fastq"

# Demultiplex
paste -d '\0' <(gunzip -kfc $a1) <(awk '{ if (NR%2+1==1) { print } else { print "" } }'  $aindex) | gzip > Princeton_read1.fq.gz  
paste -d '\0' <(gunzip -kfc $a2) <(awk '{ if (NR%2+1==1) { print } else { print "" } }'  $aindex) | gzip > Princeton_read2.fq.gz 

gunzip -kfc Princeton_read1.fq.gz | fastx_barcode_splitter.pl --eol --mismatches 1 --bcfile princeton_barcodes.txt --prefix "Princeton-P-" --suffix "-1.fq"
gunzip -kfc Princeton_read2.fq.gz | fastx_barcode_splitter.pl --eol --mismatches 1 --bcfile princeton_barcodes.txt --prefix "Princeton-P-" --suffix "-2.fq"

# Trim Reads and QC using parallel function.

function trim (){
    # Trim last 8 bp.
    seqtk trimfq -e 8 $1 | gzip > ${1%%.fq}.trimmed.fq.gz
    # Run Quality Check
    fastqc ${1%%.fq}.trimmed.fq.gz
}

export -f trim
parallel trim ::: `ls *.fq`



