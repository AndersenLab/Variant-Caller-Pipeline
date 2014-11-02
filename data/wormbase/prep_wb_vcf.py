import os, gzip

if not os.path.isfile("c_elegans.WS245.annotations.gff3.gz"):
	print "Downloading Annotation File"
	os.system("curl 'ftp://ftp.wormbase.org/pub/wormbase/releases/WS245/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS245.annotations.gff3.gz' > c_elegans.WS245.annotations.gff3.gz")

acceptable_types = ['SNP', 'point_mutation']

def parse_info(info):
	info_set = info.strip().split(";")
	ret_dict = {}
	for x in info_set:
		key, val = x.split("=")
		ret_dict[key] = val
	return ret_dict

"""
# First loop - generate strain list
with gzip.GzipFile("c_elegans.WS245.annotations.gff3.gz") as wb:
	strain_list = []
	for line in wb:
		if line.startswith("#"):
			pass
		else:
			l = line.split("\t")
			if l[2] in acceptable_types:
				info_set = parse_info(l[8])
				if "strain" in info_set:
					strains = info_set["strain"].split(",")
					for i in strains:
						if i not in strain_list and not i.startswith("VC"):
							strain_list.append(i)
							print strain_list
"""

strain_list = ['ED3040', 'ED3042', 'ED3049', 'LKC34', 'JU1401', 'MY1', 'MY6', 'JU1400', 'CB4856', 'JU312', 'JU775', 'MY16', 'AB3', 'CB4853', 'CB4854', 'ED3021', 'JU263', 'JU300', 'JU322', 'JU345', 'JU361', 'JU397', 'JU1088', 'JU1652', 'KR314', 'PX174', 'JU1171', 'MY2', 'MY14', 'CB4858', 'ED3052', 'JU258', 'JU360', 'QX1211', 'QX1216', 'OH9330', 'CX11262', 'CX11264', 'CX11268', 'CX11271', 'CX11276', 'CX11285', 'CX11294', 'CX11315', 'CX11316', 'CX11319', 'CX11321', 'DL200', 'ED3046', 'ED3048', 'JU311', 'JU313', 'JU314', 'JU316', 'JU317', 'JU318', 'JU321', 'JU393', 'JU398', 'JU401', 'JU406', 'JU622', 'JU751', 'JU778', 'JU782', 'JU799', 'JU847', 'JU848', 'JU1246', 'JU1409', 'JU1410', 'JU1411', 'JU1484', 'JU1491', 'JU533', 'AB4', 'ED3043', 'JU774', 'JU829', 'JU830', 'MY7', 'MY10', 'ED3054', 'ED3057', 'ED3072', 'JU642', 'AB2', 'CB3198', 'CB4857', 'CX11258', 'CX11259', 'CX11278', 'CX11292', 'CX11305', 'CX11307', 'CX11314', 'CX11317', 'DL226', 'DR1350', 'ED3005', 'ED3011', 'ED3077', 'JT11362', 'JT11398', 'JT11399', 'JU262', 'JU299', 'JU310', 'JU323', 'JU342', 'JU346', 'JU347', 'JU367', 'JU438', 'JU693', 'JU694', 'JU1026', 'JU1037', 'JU1039', 'JU1040', 'JU1172', 'JU1204', 'JU1206', 'JU1530', 'PB306', 'PS2025', 'RC301', 'OH7317', 'ED3017', 'OH2042', 'CB4855', 'DL238', 'EG4348', 'EG4349', 'EG4724', 'JU395', 'JU396', 'JU792', 'MY15', 'MY18', 'MY23', 'QX1218', 'QX1233', 'LSJ1', 'LSJ2', 'JU1511', 'JU1522', 'JU1581', 'JU1582', 'OH4240', 'OH4247', 'MY3', 'AB1', 'GXW1', 'OH7677', 'DM1017', 'OH8001', 'ED3073', 'ED3075', 'EG4725', 'JU531', 'JU1213', 'JU1440', 'JU1482', 'JU1580', 'JU394', 'OH8547', 'JU1218', 'JU1230', 'JU1242', 'JU1243', 'OH6071', 'IS656', 'IS659', 'OH6087', 'OH8421', 'OH8545', 'OH7116', 'OH4303', 'OH9482', 'IS657', 'OH9305', 'OH9331', 'CB4852', 'CB4932', 'EG4346', 'EG4347', 'EG4680', 'EG4689', 'EG4945', 'EG4946', 'EG4948', 'EG4951', 'EG4957', 'JU362', 'JU363', 'JU368', 'JU561', 'JU563', 'JU801', 'JU1516', 'JU1563', 'JU1896', 'PB303', 'PX178', 'PX179', 'WN2002', 'ED3019', 'ED3020', 'RB5002', 'RB5000', 'ED3010', 'JU1566', 'JU1568', 'CB4851', 'DR1344', 'RW7000', 'ED3012', 'ED3014', 'ED3015', 'ED3023', 'ED3024', 'ED3028', 'JU1395', 'WN2003', 'JU399', 'JU440', 'JU1200', 'JU1207', 'JU1212', 'JU1214', 'JU1586', 'RB5001', 'CB998', 'OH998', 'DM1154', 'CB4869', 'OH7076', 'OH7103', 'OH7075', 'OH7074', 'OH7079']


vcf_header = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=CHROMOSOME_I,length=15072423>
##contig=<ID=CHROMOSOME_II,length=15279345>
##contig=<ID=CHROMOSOME_III,length=13783700>
##contig=<ID=CHROMOSOME_IV,length=17493793>
##contig=<ID=CHROMOSOME_V,length=20924149>
##contig=<ID=CHROMOSOME_X,length=17718866>
##contig=<ID=CHROMOSOME_M,length=13794>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""

vcf = open("WS245.vcf",'w+')
vcf.write(vcf_header)
vcf.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s\n" % '\t'.join(strain_list))


def process_genotypes(strains):
	gt_set = []
	print strains
	for i in strain_list:
		if i in strains:
			gt_set.append("1/1:100")
		else:
			gt_set.append("0/0:100")
	return '\t'.join(gt_set)

c = 0
with gzip.GzipFile("c_elegans.WS245.annotations.gff3.gz") as wb:
	for line in wb:
		if c > 1000:
			break
		if line.startswith("#"):
			pass
		else:
			l = line.strip().split("\t")
			if l[2] in acceptable_types:
				info_set = parse_info(l[8])
				CHROM = "CHROMOSOME_" + l[0]
				POS = l[3]
				ID = info_set["variation"]
				TYPE = l[2]
				if "substitution" in info_set.keys() and "strain" in info_set.keys():
					variant_strains = [x for x in info_set["strain"].split(",") if not x.startswith("VC")]
					if len(variant_strains) > 0:
						GENOTYPES = process_genotypes(variant_strains)
						REF, ALT = info_set["substitution"].split("/")
						vcf_record = "{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t100\t.\t.\tGT:GQ\t{GENOTYPES}\n".format(**locals())
						vcf.write(vcf_record)
						c += 1
			else:
				pass

# Convert to vcf File
os.system("bcftools view -O z WS245.vcf > WS245.vcf.gz")
os.system("bcftools index WS245.vcf.gz")