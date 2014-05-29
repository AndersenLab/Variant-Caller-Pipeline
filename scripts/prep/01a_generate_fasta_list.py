import glob

fasta_list = glob.glob("*.fq.gz")

f = file("fastq_list.txt", "w")

for ele in zip([x for x in fasta_list if "-1.fq.gz" in x], [x for x in fasta_list if "-2.fq.gz" in x]):
	f.write("%s\t%s\n" % (ele[0], ele[1]))

f.close()