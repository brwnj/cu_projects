#!/usr/bin/env python
# encoding: utf-8
"""
Align observed CDR3 sequences to target sequences obtained outside of NGS.
"""

import os
import sys
import pandas as pd
from Bio import pairwise2
from toolshed import reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def main(queries, targets, identity, shortest):
    query_sets = {}

    for query_file in queries:
        sample_id = os.path.basename(query_file).split("_", 2)[-1].rsplit("_", 1)[0]
        queries = set([toks['CDR3-IMGT'] for toks in reader(query_file) if toks['Functionality'] == "productive" and len(toks['CDR3-IMGT']) >= shortest])
        query_sets[sample_id] = queries

    # process target sequences
    target_df = pd.read_table(targets)

    # data['patient']['peptide_source:peptide_sequence:miseq_result'] = 1
    data = {}

    for sample_id, queries in query_sets.iteritems():

        print >>sys.stderr, "processing peptides for", sample_id
        data[sample_id] = {}

        for source, series in target_df.iteritems():

            for idx, target in series.iteritems():

                if not isinstance(target, basestring): continue
                target_length = len(target)

                if target_length < shortest: continue

                for query in queries:

                    query_length = len(query)

                    longest = float(max([target_length, query_length]))

                    for a_query, a_target, score, start, stop in pairwise2.align.localms(query, target, 1, -1, -5, -1):
                        if score / longest > identity:
                            # add matched alignment
                            data[sample_id]["{source}:{known}:{miseq}".format(source=source, known=target, miseq=query)] = 1


    df = pd.DataFrame(data)
    df.index = pd.MultiIndex.from_tuples([x.split(":") for x in df.index], names=['peptide_source', 'peptide_sequence', 'miseq_result'])
    df.to_csv(sys.stdout, sep="\t", na_rep="0")

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__,
            formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('-q', '--queries', nargs="+", help="IMGT AA sequence files -- should be 5_*.txt")
    p.add_argument('-t', '--targets', help="target sequence file with peptide sequences per column, header with source id")
    p.add_argument("-i", "--identity", type=float, default=0.75, help="sequence identity threshold -- number identical divided by length of the shorter sequence")
    p.add_argument("-s", "--shortest", type=int, default=6, help="shortest allowable peptide sequence")
    args = vars(p.parse_args())

    main(**args)
