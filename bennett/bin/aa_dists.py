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


def cluster_proteins(sequences, sample, threshold, word_size):
    """
    sequences   pd.Series of AA sequences

    returns list of sequences.
    """
    # create fasta of sequences
    input_fasta = open("%s_cdhit_in.fasta" % sample, "wb")
    output_fasta = "%s_cdhit_out.fasta" % sample
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
    """

    cluster_word_size = word_size(similarity)

    print >>sys.stderr, ">>>", imgt_aa

    filename, ext = os.path.splitext(imgt_aa)

    # import into a table; rename cols
    df = pd.read_table(imgt_aa, header=0, usecols=[0,2,3,4,6,14],
                        compression="gzip" if ext == ".gz" else None,
                        names=["seq_num", "functionality", "v_gene", "j_gene", "vdj_seq", "imgt_cdr3"])

    # drop unproductive translations
    df = df[df.functionality == "productive"]
    # only take one translated germline
    df['v_gene'] = df['v_gene'].apply(lambda x: pd.Series(x.split(" ", 2)[1].split("*")[0]))
    df['j_gene'] = df['j_gene'].apply(lambda x: pd.Series(x.split(" ", 2)[1].split("*")[0]))
    # vh leader
    df['vleader'] = df['v_gene'].apply(lambda x: pd.Series(x.split("-")[0]))
    # vh exon
    # df['vexon'] = df['v_gene'].apply(lambda x: pd.Series(x.split("-")[-1]))
    # length of cdr3 sequence to group by
    df['cdr3_length'] = df.imgt_cdr3.apply(len)

    # find unique CDR3 sequences
    unique_seqs = set()
    for l, grouped_df in df.groupby('cdr3_length'):
        if l < minimum: continue

        clustered = cluster_proteins(grouped_df['imgt_cdr3'], imgt_aa, similarity, cluster_word_size)

        for seq in clustered:
            unique_seqs.add(seq)

    indexes = set()
    # filter table using unique sequences
    for seq, grouped_df in df.groupby('imgt_cdr3'):
        if not seq in unique_seqs: continue

        try:
            # most abundant V
            v = grouped_df.v_gene.value_counts().index[0]
            # most abundant J within those Vs
            j = grouped_df[grouped_df['v_gene'].isin([v])].j_gene.value_counts().index[0]
            # most abundant VDJ within the V and J group
            vdj = grouped_df[grouped_df['v_gene'].isin([v]) & grouped_df['j_gene'].isin([j])].vdj_seq.value_counts().index[0]
        except IndexError:
            # NaN in either V, J, or VDJ
            continue

        # index of row matching V, J, VDJ
        idx = grouped_df[grouped_df['v_gene'].isin([v]) & grouped_df['j_gene'].isin([j]) & grouped_df['vdj_seq'].isin([vdj])].index[0]
        indexes.add(idx)

    # no longer true because we didn't remove missing data before
    # finding unique seqs
    # assert len(unique_seqs) == len(indexes)

    # print table
    sorted_idx = sorted(indexes)
    ss = df.ix[sorted_idx]
    ss.to_csv(sys.stdout, sep="\t", cols=['vleader', 'v_gene', 'j_gene', 'vdj_seq', 'imgt_cdr3'])

    # print composition to stderr
    t = float(len(ss))
    counts = ss.vleader.value_counts()
    for k, v in counts.iteritems():
        perc = (v / t) * 100
        fields = [k, str(v), str(perc)]
        print >>sys.stderr, "\t".join(fields)
    print >>sys.stderr, "Total: %d" % t


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('imgt_aa')
    p.add_argument('-m', '--minimum', default=6, type=int, help="minimum allowable length for a productive CDR3")
    p.add_argument('-s', '--similarity', default=0.80, type=float,
        help="similarity score used to group similar CDR3 sequences: 1 \
                (sequences must be identical) to 0.4 (group almost everything \
                with the same V and J)")
    args = vars(p.parse_args())
    main(**args)
