#!/usr/local/bin/python

import os
import pandas as pd

os.chdir("../../data/ancillary/bam_sets")


os.system("ssh dec211@dhunni.biochem.northwestern.edu \"cd /lscr2/andersenlab/dec211/data/bam/ && ls *.bam  -1 | egrep -v '(AF16|HK104)'\" > 00_all_bams.txt")


# Run Sets
os.system('grep BGI1 00_all_bams.txt > 01a_BGI1_set.txt')
os.system('grep BGI2 00_all_bams.txt > 01b_BGI2_set.txt')
os.system('grep BGI3 00_all_bams.txt > 01c_BGI3_set.txt')

# BGI2 Repeats
df = pd.read_csv("00_all_bams.txt", names=['full_name'])
df['Run'], df['Library'], df['Strain'], df['hash1'], df['hash2'] = zip(*df['full_name'].map(lambda x: x.replace(".bam","").split('-')))
df = df.groupby(["Run","Library","Strain"]).filter(lambda x: x.count() == 2)

df['idx'] = df.index % 2
# Output 1st Repeat
df[df['idx'] == 0]['full_name'].to_csv('02a_BGI2_rep1.txt', index=False)
# Output 2nd Repeat
df[df['idx'] == 1]['full_name'].to_csv('02b_BGI2_rep2.txt', index=False)
# pull out dups

# Output BGI 1 & 3 for all 7 libs. 
df = pd.read_csv("00_all_bams.txt", names=['full_name'])
df['Run'], df['Library'], df['Strain'], df['hash1'], df['hash2'] = zip(*df['full_name'].map(lambda x: x.replace(".bam","").split('-')))

for lib in df['Library'].unique():
	df[df['Library']==lib]['full_name'].to_csv('03_%s.txt' % lib, index=False)

# MMP Only
mmp_strains = """AB1
AB3
CB4853
CB4854
CB4856
ED3017
ED3021
ED3040
ED3042
ED3049
ED3052
ED3057
ED3072
GXW1
JU1088
JU1171
JU1400
JU1401
JU1652
JU258
JU263
JU300
JU312
JU322
JU345
JU360
JU361
JU394
JU397
JU533
JU642
JU775
KR314
LKC34
MY1
MY14
MY16
MY2
MY6
PX174""".split('\n')

df = pd.read_csv("00_all_bams.txt", names=['full_name'])
df['Run'], df['Library'], df['Strain'], df['hash1'], df['hash2'] = zip(*df['full_name'].map(lambda x: x.replace(".bam","").split('-')))
df['full_name'][df['Strain'].isin(mmp_strains)].to_csv("04_mmp_strains.txt",index=False)


