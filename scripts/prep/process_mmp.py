# Process the vcfs

f = open("../../data/mmp/CB4856.wb225.fa.combined.vcf",'r')

for line in f.readlines()[0:500]:
	if line.startswith("##"):
		pass
	elif line.startswith("#"):
		print line
	else:
		# Action to take if line is a variant.
		variant = line.split('\t')[2].split("_")
		print variant