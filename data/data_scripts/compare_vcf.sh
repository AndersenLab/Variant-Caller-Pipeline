#!/bin/bash

# This script takes two vcf files as input and outputs a comparison of the two.

bedtools intersect -a $1 -b $2 > both.vcf