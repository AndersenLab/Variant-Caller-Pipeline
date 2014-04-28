#!/bin/bash


# Upload

cd ../raw/shipma_bam

strains="AB1
CB4854
CB4856 (CGC)
ED3017
ED3040
ED3049
ED3052
GXW1
JU1088
JU1400
JU1652
JU258
JU360
JU394
JU397
JU642
JU775
KR314
LKC34
MY1
MY16"

for r in $strains 
do
	scp `ls *$r*` dec211@dhunni:/lscr2/andersenlab/dec211/data/bam/shipma_bam/
done

PICARD="../../../tools/"



# Fix Headers
for file in `ls BGI3_RET7a_QG558.bam`
do
	f1=(${file//_/ })
	samtools view -H $file | grep '@RG'
	mv $file tmp.${file}
	java	-jar	${PICARD}AddOrReplaceReadGroups.jar	\
	I=tmp.${file}	\
	O=${file}	\
	RGID=`basename ${file}`	RGLB=${f1[1]}	\
	RGPL=Illumina RGPU=${f1[1]} \
	RGSM=${f1[2]/.bam/}	RGCN=BGI \
	VALIDATION_STRINGENCY=SILENT
	SO=coordinate
	rm tmp.${file}
done