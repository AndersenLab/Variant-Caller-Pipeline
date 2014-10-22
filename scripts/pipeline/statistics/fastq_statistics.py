#!/usr/bin/python

#=============#
# Description #
#=============#
#
# This script will submit a job that:
#
# * Generate statistics for *all* fastqs in the fq directory.

#SBATCH --job-name=fqstats

#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=16384


#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/fq

from fastq import fastq
from db_conn import * # Import sqlite db connects
import glob

def parse_fq_filename(f):
	r = dict(zip(["Run","Library","Strain","hash","pair"], f.split("-")))
	r["pair"] = r["pair"].replace(".fq.gz", "")
	return r


# Generate statistics for fastq's that do not yet have them.
files = glob.glob("*fq.gz")

fastq_set = []
# Check if stats exist within database
for i in files:
	item_count = EAV.select().where(EAV.Sub_Entity == i).count()
	if item_count == 0:
		file_info = parse_fq_filename(i)
		# Stats have not been run on the fq; run the stats.
		fq = fastq(i)
		fq.fastq_stats()
		for key,val in fq.fq1.items():
			print key, val
			save_eav(file_info["Strain"], key, val, Entity_Group = "Fastq Statistics", Sub_Entity = i[0], Tool="Awk")


