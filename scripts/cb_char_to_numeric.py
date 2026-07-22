#!/bin/python

import gzip, sys

fasta, infile, outfile = sys.argv[1:4]

bc2id = {}
with open(fasta) as fa:
    for line in fa:
        if line.startswith(">"):
            num = line[1:].strip()
        else:
            bc2id[line.strip()] = num

with open(infile) as f, gzip.open(outfile, "wt") as out:
    for line in f:
        out.write(bc2id[line.strip()] + "\n")
