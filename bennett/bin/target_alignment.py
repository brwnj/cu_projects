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


def main(queries, targets, imgt_column, identity, length):
    query_sets = {}

    for query_file in queries:
        sample_id = os.path.basename(query_file).split("_", 2)[-1].rsplit("_", 1)[0]
        query_sets[sample_id] = {}

        for toks in reader(query_file):
            if not toks['Functionality'] == "productive": continue

            # read_7:AAAAAAAT:Cm2:VH3Fr1
            vh_family = toks['Sequence ID'].rsplit(":", 1)[1]
            # only care about unique CDR3s
            query_sets[sample_id][toks[imgt_column]] = {'vdj':toks['V-D-J-REGION'], 'vh_family':vh_family}

    # process target sequences
    target_df = pd.read_table(targets)

    # data['id']['peptide_src:peptide_seq:miseq_res:vdj:vh_family'] = 1
    data = {}

    for sample_id, cdr3s in query_sets.iteritems():

        print >>sys.stderr, "processing peptides for", sample_id
        data[sample_id] = {}

        for source, series in target_df.iteritems():

            for idx, target in series.iteritems():

                if not isinstance(target, basestring): continue
                target_length = len(target)

                if target_length < shortest: continue

                for query, toks in cdr3s.iteritems():

                    query_length = len(query)

                    # only compare sequences that are similar length
                    if abs(query_length - target_length) > length: continue

                    # find shortest for identity calculation
                    shortest = float(min([target_length, query_length]))

                    # find edit distance accounting for length disparity
                    d = distance(query, target)

                    # add matched alignment
                    if ((shortest - d) / shortest) > identity:
                        k = "{source}:{known}:{miseq}:{vdj}:{vh}".format(\
                                source=source,
                                known=target,
                                miseq=query,
                                vdj=toks['vdj'],
                                vh=toks['vh_family'])
                        data[sample_id][k] = 1

    df = pd.DataFrame(data)
    df.index = pd.MultiIndex.from_tuples([x.split(":") for x in df.index], \
                names=['peptide_source', 'peptide_sequence', 'miseq_result', 'VDJ', 'VH'])
    df.to_csv(sys.stdout, sep="\t", na_rep="0")


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__,
            formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('-q', '--queries', nargs="+", help="IMGT AA sequence files -- should be 5_*.txt")
    p.add_argument('-t', '--targets', help="target sequence file with peptide sequences per column, header with source id")
    p.add_argument("-c", "--imgt-column", default="CDR3-IMGT", help="column name for sequences")
    p.add_argument("-i", "--identity", type=float, default=0.75, help="sequence identity threshold -- number identical divided by length of the shorter sequence")
    p.add_argument("-l", "--length", type=int, default=2, help="allowable length difference between target and query sequence")
    args = vars(p.parse_args())

    main(**args)
