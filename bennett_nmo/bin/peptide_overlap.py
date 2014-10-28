#!/usr/bin/env python
# encoding: utf-8
"""
Find peptide overlaps between CSF and Blood CDR3 sequences.
"""
import os
import sys
import argparse
import pandas as pd
import subprocess as sp
from toolshed import nopen
from itertools import groupby

def readfa(fa):
    with nopen(fa) as fh:
        for header, group in groupby(fh, lambda line: line[0] == '>'):
            if header:
                line = group.next()
                name = line[1:].strip()
            else:
                seq = ''.join(line.strip() for line in group)
                yield name, seq

def csf_to_set(filename):
    s = set()
    for name, seq in readfa(filename):
        s.add(seq)
    return s

def run_cd_hit_2d(fasta_a, fasta_b, threshold):
    """
    cd-hit-2d output will contain sequences of novel sequences of fasta_b
    """
    output_fasta = "cdhit2d_out.fasta"
    # run fasta through cd-hit
    cmd = "cd-hit-2d -i {input1} -i2 {input2} -o {output} -c {threshold} -d 0".format(input1=fasta_a, input2=fasta_b, output=output_fasta, threshold=threshold)
    p = sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    # parse output for remaining sequences
    novel = [seq for name, seq in readfa(output_fasta)]
    # cleanup files from cd-hit call
    os.remove(output_fasta)
    os.remove("{output}.clstr".format(output=output_fasta))
    # return set of non-overlapping peptides
    return set(novel)

def main(csf, blood, threshold):
    # convert csf into set
    csf_set = csf_to_set(csf)
    # run cd-hit-2d for each blood file
    peptides = {}
    for fasta in blood:
        sample = fasta.split(".", 1)[0]
        peptides[sample] = {}
        novel_peptides = run_cd_hit_2d(fasta, csf, threshold)
        overlapping_peptides = csf_set - novel_peptides
        for peptide in overlapping_peptides:
            peptides[sample][peptide] = True
    df = pd.DataFrame(peptides)
    df.to_csv(sys.stdout, sep="\t", na_rep=False)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('csf', help="csf peptide fasta")
    p.add_argument('blood', nargs="+", help="blood peptide fasta(s)")
    p.add_argument('-t', '--identity-threshold', type=float, default=0.90,
            help="Sequence Identity threshold: the number of identical amino \
            acids in alignment divided by the full length of the shorter \
            sequence.")
    args = p.parse_args()
    main(args.csf, args.blood, args.identity_threshold)
