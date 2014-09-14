#!bin/python

import gzip
import os

with gzip.open('c_elegans.WS220.genomic_masked.fa.gz', 'rb') as f:
    file_content = f.read()
    print file_content