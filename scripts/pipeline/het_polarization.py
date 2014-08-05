
'''
Het Polarization Script
usage:
bcftools view <filename> | head -n 500 | python het_polarization.py  | bcftools view -O b > 00_ALL_hets_polarized
'''

import sys

log = open("het_polarization_log.txt",'w')

PL_diff = 20

def main():
    for l in sys.stdin.xreadlines():
        het_set = []
        if l.startswith("#CHROM"):
            # Get Sample information and count
            samples = l.split("\t")[9:]
            log.write("CHROM\tPOS\t" + '\t'.join(samples))
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
                    PL_set = [int(i) for i in v.split(":")[PL].split(",")]
                    if (PL_set[0] - PL_set[2]) > PL_diff:
                        l[k+9] = v.replace("0/1","0/0")
                    elif (PL_set[2] - PL_set[0]) > PL_diff:
                        l[k+9] = v.replace("0/1","1/1")
                    else:
                        pass
                    het_set.append((PL_set[2] - PL_set[0]))
                else:
                    het_set.append(0)
            sys.stdout.write("\t".join(l))
            log.write("\t".join(map(str,het_set)) + "\n")
        else:
            het_set.append(l[0:l.index("\t",l.index("\t")+1)] + '\t'.join((1+len(samples))*['0']))
            sys.stdout.write(l)
            log.write("\t".join(map(str,het_set)) + "\n")

if __name__ == '__main__':
    main()