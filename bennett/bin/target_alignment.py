#!/usr/bin/env python
# encoding: utf-8
"""
Align observed CDR3 sequences to target sequences obtained outside of NGS.
"""

import os
import sys
import editdist
import pandas as pd
from toolshed import reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def distance(a, b):
    """
    Find best edit distance between two strings of potentially uneven length.

    >>> import editdist
    >>> distance("abc", "abc")
    0
    >>> distance("abc", "abcdef")
    0
    >>> distance("abbb", "abcdef")
    2
    """
    # first position longer
    ab = sorted([a, b], key=len, reverse=True)
    assert all(isinstance(s, basestring) for s in ab)

    length_longer, length_shorter = len(ab[0]), len(ab[1])

    if length_longer == length_shorter:
        d = editdist.distance(ab[0], ab[1])

    else:
        dists = set()
        difference = length_longer - length_shorter

        for i in xrange(0, difference + 1):
            d = editdist.distance(ab[0][i:i + length_shorter], ab[1])
            # perfectly matching substring
            if d == 0:
                break
            dists.add(d)

        if d != 0:
            d = min(dists)

    return d


def match_scores(queries, targets, length):
    # pull out all unique CDR3 sequences
    blood_cdr3s = set()
    for query_file in queries:
        for t in reader(query_file):
            if t['Functionality'] == 'productive':
                blood_cdr3s.add(t['CDR3-IMGT'])

    # calculate edit distances for all unique CDR3s; trie would be ideal here
    target_df = pd.read_table(targets)
    alignments = {}
    seq_to_source = {}
    for query in blood_cdr3s:
        query_length = len(query)
        alignments[query] = {}

        for source, series in target_df.iteritems():
            # not sure if this actually works as intended
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

    return alignments, seq_to_source


def process_queries(queries, alignments, seq_to_source, identity):
    """

    queries       - AA IMGT result files
    alignments    - dict of dicts; blood CDR3 as key, CSF CDR3 as key, match score as value
    seq_to_source - dict; sequence as key, peptide source as value
    identity      - threshold for sequence similarity
    """
    for query_file in queries:
        sample = query_file.split("sequences_")[1].rsplit("_", 1)[0]
        print >>sys.stderr, "processing", sample
        outfile = open(sample + "_aligned.txt", 'w')

        # import into a table
        df = pd.read_table(query_file, index_col=0, header=0, usecols=[0,1,2,3,4,6,14],
                            names=["seq_num", "seq_id", "functionality", "v_gene",
                                    "j_gene", "vdj_seq", "imgt_cdr3"])

        # drop unproductive translations
        df = df[df.functionality == "productive"]
        # only take one translated germline
        df['v_gene'] = df['v_gene'].apply(lambda x: pd.Series(x.split(" ", 2)[1]))
        df['j_gene'] = df['j_gene'].apply(lambda x: pd.Series(x.split(" ", 2)[1]))
        # primer germline
        df['primer_germline'] = df['seq_id'].apply(lambda x: pd.Series(x.rsplit(":", 1)[1]))
        df['read_id'] = df['seq_id'].apply(lambda x: pd.Series(x.split(":", 1)[0]))
        df['umi'] = df['seq_id'].apply(lambda x: pd.Series(x.split(":")[1]))

        # print the sequences
        header = ["read_id", "umi", "number_collapsed", "csf_cdr3", "imgt_cdr3",
                    "match_score", "csf_source", "primer_germline", "v_gene", "vdj_seq"]
        print >>outfile, "\t".join(header)

        # groupby imgt_cdr3, primer_germline, v_gene, vdj_seq
        for (imgt_cdr3, primer_germline, v_gene, vdj_seq), grouped_df in \
                df.groupby(['imgt_cdr3', 'primer_germline', 'v_gene', 'vdj_seq']):

            # compile list of read_ids and umis
            reads = list(grouped_df['read_id'])
            umis = list(grouped_df['umi'])
            collapsed = len(reads)
            if collapsed > 10:
                reads = reads[:10]
                umis = umis[:10]

            # look up all alignments, print where above threshold
            for csf_cdr3, match_score in alignments[imgt_cdr3].iteritems():
                if match_score > identity:
                    fields = [",".join(reads), ",".join(umis), collapsed, csf_cdr3,
                                imgt_cdr3, match_score, seq_to_source[csf_cdr3],
                                primer_germline, v_gene, vdj_seq]
                    print >>outfile, "\t".join([str(f) for f in fields])

        outfile.close()


def main(queries, targets, identity, length):
    """

    queries  - aa-sequence result files from IMGT
    targets  - aa-sequence results from markus
    identity - match score threshold
    length   - length allowance when calculating match_score
    """

    print >>sys.stderr, "calculating match scores for all observed peptides sequences"
    alignments, seq_to_source = match_scores(queries, targets, length)
    process_queries(queries, alignments, seq_to_source, identity)


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__,
            formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('-q', '--queries', nargs="+", help="IMGT AA sequence files -- should be 5_*.txt")
    p.add_argument('-t', '--targets', help="target sequence file with peptide sequences per column, header with source id")
    p.add_argument('-i', '--identity', type=float, default=0.75, help="sequence identity threshold -- number identical divided by length of the shorter sequence")
    p.add_argument('-l', '--length', type=int, default=2, help="allowable length difference between target and query sequence -- -1 to disable")
    args = vars(p.parse_args())

    main(**args)
