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


def cluster_proteins(sequences, threshold, word_size):
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
            "-g 1 -t 1 -l 4 -M 4000 -n {word_size}").format(input=input_fasta.name,
                                                  output=output_fasta,
                                                  threshold=threshold,
                                                  word_size=word_size)

    p = sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    # parse output for remaining sequences
    clustered_proteins = [seq for name, seq in readfa(output_fasta)]
    # cleanup files from cd-hit call
    rm([input_fasta.name, output_fasta, output_fasta + ".clstr"])
    # return list of remaining sequences
    return clustered_proteins


def word_size(n):
    assert 1.0 >= n >= 0.4
    if n > 0.7:
        ws = 5
    elif n > 0.6:
        ws = 4
    elif n > 0.5:
        ws = 3
    else:
        ws = 2
    return ws


def main(imgt_aa, minimum, similarity):
    """
    + same v and j region, unlikely to have different CDR3 sequence
    + need to have another column in the metadata counting those unique sequences in this manner
    + don't care about the actual CDR3 sequence at this point
    + need population dist for these sequences as well
    """
    cluster_word_size = word_size(similarity)

    imgt_aa = "/Users/brownj/projects/bennett/data/common/ON10_03B_memory/5_AA-sequences_ON10_03B_memory_170114.txt"
    minimum = 6
    similarity = 0.8


    sample = imgt_aa.split("sequences_")[1].rsplit("_", 1)[0]
    print >>sys.stderr, "processing", sample
    # outfile = open(sample + "_aligned.txt", 'w')

    # import into a table; rename cols
    df = pd.read_table(imgt_aa, index_col=0, header=0, usecols=[0,2,3,4,14],
                        names=["seq_num", "functionality", "v_gene", "j_gene", "imgt_cdr3"])

    # drop unproductive translations
    df = df[df.functionality == "productive"]
    # only take one translated germline
    df['v_gene'] = df['v_gene'].apply(lambda x: pd.Series(x.split(" ", 2)[1].split("*")[0]))
    df['j_gene'] = df['j_gene'].apply(lambda x: pd.Series(x.split(" ", 2)[1].split("*")[0].strip("IGH")))
    # vh leader
    df['vleader'] = df['v_gene'].apply(lambda x: pd.Series(x.strip("IGH").split("-")[0]))
    # vh exon
    df['vexon'] = df['v_gene'].apply(lambda x: pd.Series(x.split("-")[-1]))
    # length of cdr3 sequence to group by
    df['cdr3_length'] = df.imgt_cdr3.apply(len)

    # meeting notes
    # similarity should be very high, like .9
    # then put into groups for quantification
    # output the table of of cdr3 and v_gene name and j_gene name

    # find unique CDR3 sequences
    unique_seqs = set()
    for l, grouped_df in df.groupby('cdr3_length'):
        if l < minimum: continue

        clustered = cluster_proteins(grouped_df['imgt_cdr3'], similarity, cluster_word_size)

        for seq in clustered:
            unique_seqs.add(seq)

    # filter table using unique sequences
    for seq, grouped_df in df.groupby('imgt_cdr3'):
        if not seq in unique_seqs: continue
        # find most abundant v:j in grouped_df

        # grouped_df is full of the same CDR3 sequence
        # within grouped_df find most abundant entry

# for s,t in df.groupby('imgt_cdr3'):
#     print s
#     v = t.v_gene.value_counts().index[0]
#     j = t.j_gene.value_counts().index[0]
#     idx = t[t['v_gene'].isin([v])].index[0]
#     print t.ix[idx]
#     assert t.ix[idx].j_gene == j
#     break

        # most abundant v-gene; want whole entry
        v = grouped_df.v_gene.value_counts().index[0]
        j = grouped_df.j_gene.value_counts().index[0]


    df = pd.DataFrame(counts)
    df.to_pickle("something")


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('imgt_aa')
    p.add_argument('-m', '--minimum', default=6, type=int, help="minimum allowable length for a productive CDR3")
    p.add_argument('-s', '--similarity', default=0.80, type=float,
        help="similarity score used to group similar CDR3 sequences: 1 \
                (sequences must be identical) to 0 (just group everything \
                with the same V and J)")
    args = vars(p.parse_args())
    main(**args)
