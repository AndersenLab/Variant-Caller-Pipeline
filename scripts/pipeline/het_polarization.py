'''
Heterozygote Polarization Script
usage:
bcftools view <filename> | python het_polarization.py  | bcftools view -O b > <filename.het.polarized.bcf>
'''
 
import sys
 
log = open("het_polarization_log.txt",'w')
 
PL_diff = 20
 
het_prior = 0.000001
 
def phred2p(phred):
    return 10**(phred/-10.0)
 
def main():
    for l in sys.stdin.xreadlines():
        het_set = []
        if l.startswith("#CHROM"):
            # Get Sample information and count
            samples = l.replace("\n","").split("\t")[9:]
            log = open("het_polarization_log.%s.txt" % '.'.join(samples),'w')
            log.write("CHROM\tPOS\t" + '\t'.join(samples) + "\n")
            sys.stdout.write(l)
        elif l.startswith("#"):
            # Pass comment lines.
            sys.stdout.write(l)
        elif l.find("0/1") != -1:
        # Check and see if there are any hets
            l = l.split("\t")
            het_set = [l[0],l[1]]
            PL = l[8].split(":").index("PL")
            for k,v in enumerate(l[9:]):
                if v.startswith("0/1"):
                    PL_set = [phred2p(int(i)) for i in v.split(":")[PL].split(",")]
                    cond_prob = [PL_set[0]*((1-het_prior)/2), PL_set[1]*het_prior, PL_set[2]*(1-het_prior)/2]
                    cond_prob = [cond_prob[0]/sum(cond_prob), cond_prob[1]/sum(cond_prob), cond_prob[2]/sum(cond_prob)]
                    if (max(cond_prob) == cond_prob[0]):
                        l[k+9] = v.replace("0/1","0/0")
                        het_set.append("AA - %s" % cond_prob[0])
                    elif (max(cond_prob) == cond_prob[2]):
                        l[k+9] = v.replace("0/1","1/1")
                        het_set.append("BB - %s" % cond_prob[2])
                    else:
                        het_set.append("AB - %s" % cond_prob[1])
                else:
                    het_set.append(0)
            sys.stdout.write("\t".join(l))
            log.write("\t".join(map(str,het_set)) + "\n")
        else:
            het_set.append(l[0:l.index("\t",l.index("\t")+1)] + '\t'.join((1+len(samples))*['0']))
            sys.stdout.write(l)
 
if __name__ == '__main__':
    main()


