#!bin/python

import gzip
import io

with io.TextIOWrapper(io.BufferedReader(gzip.open('c_elegans.WS220.genomic_masked.fa.gz'))) as f:
    l = 0 # current line
    n = 0 # Current character
    state = 0 # 0=Out of gap; 1=In Gap
    for line in f:
        line = line.replace("\n","")
        if line.startswith(">"):
            chr = line
        else:
            if l < 100000:
                for char in line:
                    n += 1
                    if state == 0 and char == "N":
                        state = 1
                        start = n
                    elif state == 1 and char != "N":
                        state = 0
                        end = n-1
                        print chr, "|" ,start, end
                    else:
                        pass
                    #print char, n, state
                l += 1
            #print start, end
