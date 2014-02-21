#!/usr/bin/env python
# encoding: utf-8
"""
Align observed CDR3 sequences to target sequences obtained outside of NGS.
"""

import os
import sys
import pandas as pd
import editdist as ed
from toolshed import reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def distance(a, b):
    """
    find best edit distance between two strings of potentially uneven length.
    """
    la, lb = len(a), len(b)
    assert isinstance(a, basestring), isinstance(b, basestring)
    if la < lb:
        return distance(b, a)
    if la == lb:
        return ed.distance(a, b)
    else:
        dists = []
        for i in xrange(0, la-lb+1):
            dists.append(ed.distance(a[i:i+lb], b))
        return min(dists)


def main(queries, targets, identity, length):
    """
    read in all aa sequences from one sample into table
    run alignments
    print output as separate files
    """

    query = queries[0]

    # import into a table
    df = pd.read_table(query, index_col=0, header=0, usecols=[0,1,2,3,4,6,14],
                        names=["seq_num", "seq_id", "functionality", "v_gene",
                                "j_gene", "vdj_seq", "imgt_cdr3"])

    print >>sys.stderr, "munging query data for", query
    # drop unproductive translations
    df = df[df.functionality == "productive"]
    # only take one translated germline
    df['v_gene'] = df['v_gene'].apply(lambda x: pd.Series(x.split(" ", 2)[1]))
    df['j_gene'] = df['j_gene'].apply(lambda x: pd.Series(x.split(" ", 2)[1]))
    # primer germline
    df['primer_germline'] = df['seq_id'].apply(lambda x: pd.Series(x.rsplit(":", 1)[1]))
    df['read_id'] = df['seq_id'].apply(lambda x: pd.Series(x.split(":", 1)[0]))
    df['umi'] = df['seq_id'].apply(lambda x: pd.Series(x.split(":")[1]))

    # calculate edit distances for all unique CDR3s
    print >>sys.stderr, "calculating match scores"
    target_df = pd.read_table(targets)
    alignments = {}
    seq_to_source = {}
    for query, grouped_df in df.groupby('imgt_cdr3'):
        query_length = len(query)
        alignments[query] = {}
        for source, series in target_df.iteritems():
            series = series.drop_duplicates()
            for idx, target in series.iteritems():
                # columns contain uneven number of rows
                if not isinstance(target, basestring): continue

                # handle duplicates but fail on identical from separate sources
                if target in seq_to_source and source != seq_to_source[target]:
                    print >>sys.stderr, "ambigous target", target, source
                    sys.exit(1)

                seq_to_source[target] = source

                target_length = len(target)

                # length too long
                if length != -1 and abs(query_length - target_length) > length:
                    alignments[query][target] = 0
                    continue

                # length in range
                shorter = float(min([target_length, query_length]))
                d = distance(query, target)
                match_score = (shorter - d) / shorter
                alignments[query][target] = match_score

    # print the sequences
    header = ["read_id", "umi", "csf_cdr3", "imgt_cdr3", "match_score",
                "csf_source", "primer_germline", "v_gene", "vdj_seq"]

    # groupby imgt_cdr3, primer_germline, v_gene, vdj_seq
    for (imgt_cdr3, primer_germline, v_gene, vdj_seq), (dframe) in \
            df.groupby(['imgt_cdr3', 'primer_germline', 'v_gene', 'vdj_seq']):

        # compile list of read_ids and umis
        reads = list(dframe['read_id'])
        umis = list(dframe['umi'])

        # look up all alignments, print where above threshold
        for csf_cdr3, match_score in alignments[imgt_cdr3].iteritems():
            if match_score > identity:
                fields = [",".join(reads), ",".join(umis), csf_cdr3, imgt_cdr3,
                            match_score, seq_to_source[csf_cdr3],
                            primer_germline, v_gene, vdj_seq]
                print "\t".join([str(f) for f in fields])


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__,
            formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('-q', '--queries', nargs="+", help="IMGT AA sequence files -- should be 5_*.txt")
    p.add_argument('-t', '--targets', help="target sequence file with peptide sequences per column, header with source id")
    p.add_argument('-i', '--identity', type=float, default=0.75, help="sequence identity threshold -- number identical divided by length of the shorter sequence")
    p.add_argument('-l', '--length', type=int, default=2, help="allowable length difference between target and query sequence -- -1 to disable")
    args = vars(p.parse_args())

    main(**args)
