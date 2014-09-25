'''
Heterozygote Polarization Script
usage:
bcftools view <filename> | python het_polarization.py  | bcftools view -O b > <filename.het.polarized.bcf>
'''
 
import sys

het_prior = 0.000001
 
def phred2p(phred):
    return 10**(phred/-10.0)
 
def main():
    polarize_alt = 0
    polarize_ref = 0
    remain_het = 0
    info_added = False
    for l in sys.stdin.xreadlines():
        if l.startswith("#CHROM"):
            # Get Sample information and count
            samples = l.replace("\n","").split("\t")[9:]
            shortlog = open("het_polarization_log.shortlog.%s.txt" % '.'.join(samples),'w')
            sys.stdout.write(l)
        elif l.startswith("#"):
            # Add Info line for het polarization flag
            if l.startswith("##INFO") and info_added == False:
                info_added = True
                sys.stdout.write("##INFO=<ID=HetPol,Type=String,Description=\"Flag used to mark whether a variant was polarized\">\n")
            # Pass comment lines.
            sys.stdout.write(l)
        # Ignore indels, and see if there are any hets in the line before proceeding.
        elif l.find("0/1") != -1 and l.find("INDEL") == -1:
        # Check and see if there are any hets
            l = l.split("\t")
            PL = l[8].split(":").index("PL")
            for k,v in enumerate(l[9:]):
                # Exclude any line 
                if v.startswith("0/1"):
                    PL_set = [phred2p(int(i)) for i in v.split(":")[PL].split(",")]
                    cond_prob = [PL_set[0]*((1-het_prior)/2), PL_set[1]*het_prior, PL_set[2]*(1-het_prior)/2]
                    cond_prob = [cond_prob[0]/sum(cond_prob), cond_prob[1]/sum(cond_prob), cond_prob[2]/sum(cond_prob)]
                    if (max(cond_prob) == cond_prob[0]):
                        l[k+9] = v.replace("0/1","0/0")
                        l[7] = l[7] + ";HetPol=AA"
                        polarize_ref += 1
                    elif (max(cond_prob) == cond_prob[2]):
                        l[k+9] = v.replace("0/1","1/1")
                        l[7] = l[7] + ";HetPol=BB"
                        polarize_alt += 1
                    else:
                        remain_het += 1
            sys.stdout.write("\t".join(l))
        else:
            sys.stdout.write(l)
    shortlog.write("Het Call Polarized to Alt:%s,Het Call Polarized to Ref:%s,Het Call No Change:%s" % (polarize_alt, polarize_ref, remain_het))

if __name__ == '__main__':
    main()


