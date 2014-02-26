#!/usr/bin/env python
# coding=utf-8
"""
Grouping by V and J, compile distributions of CDR3 sequences.
"""

import os
import sys
import itertools
import pandas as pd
import subprocess as sp
from toolshed import nopen
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def readfa(fa):
    with nopen(fa) as fh:
        for h, g in itertools.groupby(fh, lambda l: l[0] == '>'):
            if h:
                l = g.next()
                name = l[1:].strip()
            else:
                seq = ''.join(l.strip() for l in g)
                yield name, seq

def rm(file):
    def _remove(f):
        try:
            os.remove(f)
        except OSError:
            pass

    if isinstance(file, list):
        for f in file:
            _remove(f)
    else:
        _remove(f)


def cluster_proteins(sequences, threshold):
    """
    sequences   pd.Series of AA sequences

    returns list of sequences.
    """
    # create fasta of sequences
    input_fasta = open("cdhit_in.fasta", "wb")
    output_fasta = "cdhit_out.fasta"
    # generate fasta of sequences
    for i, seq in enumerate(sequences):
        print >>input_fasta, ">seq_{id}\n{sequence}".format(id=i, sequence=seq)
    input_fasta.close()
    # run fasta through cd-hit
    cmd = ("cd-hit -i {input} -o {output} -c {threshold} "
            "-g 1 -t 1 -l 4 -M 4000 -n 2").format(input=input_fasta.name,
                                                  output=output_fasta,
                                                  threshold=threshold)

    p = sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    # parse output for remaining sequences
    clustered_proteins = [seq for name, seq in readfa(output_fasta)]
    # cleanup files from cd-hit call
    rm([input_fasta.name, output_fasta, output_fasta + ".clstr"])
    # return list of remaining sequences
    return clustered_proteins


def main(imgt_aa, similarity):
    """
    + same v and j region, unlikely to have different CDR3 sequence
    + need to have another column in the metadata counting those unique sequences in this manner
    + don't care about the actual CDR3 sequence at this point
    + need population dist for these sequences as well
    """
    sample = imgt_aa.split("sequences_")[1].rsplit("_", 1)[0]
    print >>sys.stderr, "processing", sample
    # outfile = open(sample + "_aligned.txt", 'w')

    # import into a table; rename cols
    df = pd.read_table(imgt_aa, index_col=0, header=0, usecols=[0,2,3,4,14],
                        names=["seq_num", "functionality", "v_gene", "j_gene", "imgt_cdr3"])

    # drop unproductive translations
    df = df[df.functionality == "productive"]
    # only take one translated germline
    df['v_gene'] = df['v_gene'].apply(lambda x: pd.Series(x.split(" ", 2)[1]))
    df['j_gene'] = df['j_gene'].apply(lambda x: pd.Series(x.split(" ", 2)[1]))

    total_cdr3_count = 0

    # tracking counts; thinking about pulling in all data across samples
    counts = {}
    counts[sample] = {}

    # groupby v_gene and j_gene
    for (v_gene, j_gene), grouped_df in \
            df.groupby(['v_gene', 'j_gene']):

        # need to count all unique CDR3s
        # need to count, per v-j, unique CDR3s
        sequences = cluster_proteins(grouped_df['imgt_cdr3'], similarity)
        counts[sample][v_gene + ":" + j_gene] = len(sequences)

    df = pd.DataFrame(counts)
    df.to_pickle("something")


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('imgt_aa')
    p.add_argument('-s', '--similarity', default=0.60, type=float,
        help="similarity score used to group similar CDR3 sequences: 1 \
                (sequences must be identical) to 0 (just group everything \
                with the same V and J)")
    args = vars(p.parse_args())
    main(**args)
