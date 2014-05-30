# Create a file with the princeton barcodes
echo """AF16	AGGTACC
CB4856_CGC	ACCAACT
DL238	AGTCAGA
HK104	TGACGTC
JU258	AATTCAA
MY23	CTGCGAC
N2_CGC	TCGCAGG
QX1430	CATCCGG""" >> princeton_barcodes.txt

fastx_barcode_splitter.pl princeton.fq --eol --bcfile princeton_barcodes.txt "--prefix Princeton-P-"


# Quality Reports
fastx_quality_stats -i BC54.fq -o bc54_stats.txt
fastq_quality_boxplot_graph.sh -i bc54_stats.txt -o bc54_quality.png -t "My Library"
fastx_nucleotide_distribution_graph.sh -i bc54_stats.txt -o bc54_nuc.png -t "My Library"