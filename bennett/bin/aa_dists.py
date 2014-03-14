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


def readfq(fq):
    with nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            rd = [x for x in islice(fqclean, 4)]
            if not rd: raise StopIteration
            assert all(rd) and len(rd) == 4
            yield rd[0][1:], rd[1], rd[3]


def fq_to_dict(fq):
    d = {}
    for name, seq, qual in readfq(fq):
        d[name] = seq
    return d


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


def main(imgt_aa, joined_fastq, minimum, similarity):
    """
    In naive populations, we expect 40% V3, 20% V1 and V4, low percentages for
    V2 and V5.
    """

    cluster_word_size = word_size(similarity)

    print >>sys.stderr, ">>>", imgt_aa

    filename, ext = os.path.splitext(imgt_aa)
    nt_seqs = fq_to_dict(joined_fastq)

    # import into a table; rename cols
    col_idxs = [0,1,2,3,4,6,14]
    col_names = ["seq_num", "seq_id", "functionality", "v_gene", "j_gene", "vdj_seq", "imgt_cdr3"]
    df = pd.read_table(imgt_aa, header=0, usecols=col_idxs,
                        compression="gzip" if ext == ".gz" else None,
                        names=col_names)

    # drop unproductive translations
    df = df[df.functionality == "productive"]
    # drop missing
    df.dropna(subset=['vdj_seq', 'v_gene', 'j_gene'], how='any', inplace=True)
    # only take one translated germline
    df['v_gene'] = df['v_gene'].str.split().str[1].str.split("*").str[0]
    df['j_gene'] = df['j_gene'].str.split().str[1].str.split("*").str[0]
    df['vleader'] = df['v_gene'].str.split("-").str[0]
    df['immunoglobulin'] = df['seq_id'].str.split(":").str[2].str.split("C").str[1].str.split("2").str[0]

    # length of cdr3 sequence to group by
    df['cdr3_length'] = df.imgt_cdr3.apply(len)
    # replace stars within VDJ sequence
    df['vdj_seq'] = df['vdj_seq'].str.replace("*","")

    # find unique CDR3 sequences
    unique_seqs = set()
    # added v and j gene to groupby may or may not increase len of unique seqs
    for (l, vgene, jgene), grouped_df in df.groupby(['cdr3_length', 'v_gene', 'j_gene']):
        if l < minimum: continue

        clustered = cluster_proteins(grouped_df['imgt_cdr3'], imgt_aa, similarity, cluster_word_size)

        for seq in clustered:
            unique_seqs.add(seq)

    indexes = set()
    # filter table using unique sequences
    # TODO: will need plots and tables for that as well
    for (cdr3, v_gene, j_gene), grouped_df in df.groupby(['imgt_cdr3', 'v_gene', 'j_gene']):
        if not cdr3 in unique_seqs: continue

        # prior to adding v and j in the groupby

        # most abundant V
        # v = grouped_df.v_gene.value_counts().index[0]
        # most abundant J within those Vs
        # j = grouped_df[grouped_df['v_gene'].isin([v])].j_gene.value_counts().index[0]
        # most abundant VDJ within the V and J group
        # vdj = grouped_df[grouped_df['v_gene'].isin([v]) & grouped_df['j_gene'].isin([j])].vdj_seq.value_counts().index[0]

        # index of row matching V, J, VDJ
        # idx = grouped_df[grouped_df['v_gene'].isin([v]) & grouped_df['j_gene'].isin([j]) & grouped_df['vdj_seq'].isin([vdj])].index[0]

        # most abundant VDJ sequence
        vdj = grouped_df.vdj_seq.value_counts().index[0]
        idx = grouped_df[grouped_df['vdj_seq'].isin([vdj])].index[0]

        indexes.add(idx)

    # indexes will likely be longer than unique_seqs due to grouping

    # print table
    # will need
    sorted_idx = sorted(indexes)
    ss = df.ix[sorted_idx]
    # add nucleotide sequences onto this dataframe
    ss['nt_seq'] = nt_seqs[ss['seq_id'].str]

    ss.to_csv(sys.stdout, sep="\t", cols=['immunoglobulin', 'vleader', 'v_gene', 'j_gene', 'vdj_seq', 'imgt_cdr3', 'nt_seq'])

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
    p.add_argument('joined_fastq')
    p.add_argument('-m', '--minimum', default=6, type=int, help="minimum allowable length for a productive CDR3")
    p.add_argument('-s', '--similarity', default=0.80, type=float,
        help="similarity score used to group similar CDR3 sequences: 1 \
                (sequences must be identical) to 0.4 (group almost everything \
                with the same V and J)")
    args = vars(p.parse_args())
    main(**args)
