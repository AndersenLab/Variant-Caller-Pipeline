#!/usr/bin/python

#=============#
# Description #
#=============#
#
# This script will submit a job that:
#
# * Runs reports
#



#SBATCH --job-name=bcfcomp

#SBATCH --output=../log/%j.txt
#SBATCH --error=../log/%j.out

#SBATCH --partition=compute

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem=16384

#SBATCH --mail-user=dec@u.northwestern.edu
#SBATCH --workdir=/lscr2/andersenlab/dec211/data/vcf



from datetime import datetime
startTime = datetime.now()
import os, sys

job_id = os.system("echo $SLURM_JOB_ID")

print sys.argv
print "Job Id: %s" % job_id
print
print "Command: %s" % (sys.argv)
print

# Determine type of qc being performed. 1 file, 2 files (compare) or comparison among many files (multiple)

if len(sys.argv[1:]) == 1:
	# Produce stats for single file

elif len(sys.argv[1:]) == 2:
	# Genotype check
	p=file("../reports/plot_commands.txt","a+")
	p.write("scp dec211@dhunni.biochem.northwestern.edu:/lscr2/andersenlab/dec211/data/reports/%s.stats.txt %s.stats.txt\n" % (outfile, outfile))
	p.write("scp dec211@dhunni.biochem.northwestern.edu:/lscr2/andersenlab/dec211/data/vcf/%s %s\n" % (sys.argv[2],sys.argv[2]))
	p.write("scp dec211@dhunni.biochem.northwestern.edu:/lscr2/andersenlab/dec211/data/vcf/%s %s\n" % (sys.argv[3],sys.argv[3]))
	p.write("scp dec211@dhunni.biochem.northwestern.edu:/lscr2/andersenlab/dec211/data/vcf/%s.tbi %s.tbi\n" % (sys.argv[2],sys.argv[2]))
	p.write("scp dec211@dhunni.biochem.northwestern.edu:/lscr2/andersenlab/dec211/data/vcf/%s.tbi %s.tbi\n" % (sys.argv[3],sys.argv[3]))
	p.write('bcftools gtcheck -H -p %s/ -g %s \n' % (outfile, " ".join(sys.argv[2:])))
	p.write("plot-vcfstats --main-title %s -t %s -t %s -r -p %s/ %s.stats.txt\n" % (outfile, sys.argv[2].replace(".txt","").replace(".vcf.gz",""), sys.argv[3].replace(".txt","").replace(".vcf.gz",""),outfile, outfile))
	p.write("\n\n")
else:
	print "multi"

outfile = sys.argv[1]

os.system('bcftools stats -s - %s > ../reports/%s.stats.txt' % (" ".join(sys.argv[2:]), outfile))



# Print total time
print("Time to complete: %s" % str(datetime.now()-startTime))

