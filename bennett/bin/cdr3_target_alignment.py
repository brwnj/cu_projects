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
        query_sets[sample_id] = {}

        for toks in reader(query_file):
            if not toks['Functionality'] == "productive": continue
            if not len(toks['CDR3-IMGT']) >= shortest: continue

            # read_7:AAAAAAAT:Cm2:VH3Fr1
            vh_family = toks['Sequence ID'].rsplit(":", 1)[1]
            # only care about unique CDR3s
            query_sets[sample_id][toks['CDR3-IMGT']] = {'vdj':toks['V-D-J-REGION'], 'vh_family':vh_family}

    # process target sequences
    target_df = pd.read_table(targets)

    # data['patient']['peptide_source:peptide_sequence:miseq_result:alignment:vdj:vh_family'] = 1
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

                    shortest = float(min([target_length, query_length]))

                    for a_query, a_target, score, start, stop in pairwise2.align.localms(query, target, 1, -1, -5, -1):
                        if score / shortest > identity:
                            # add matched alignment
                            # when new patient shows up, will likely not overlap
                            # do to all these things being in the key
                            k = "{source}:{known}:{miseq}:{alignment}:{vdj}:{vh}".format(\
                                    source=source,
                                    known=target,
                                    miseq=query,
                                    alignment=a_query,
                                    vdj=toks['vdj'],
                                    vh=toks['vh_family'])
                            data[sample_id][k] = 1

    df = pd.DataFrame(data)
    df.index = pd.MultiIndex.from_tuples([x.split(":") for x in df.index], \
                names=['peptide_source', 'peptide_sequence', 'miseq_result', 'alignment', 'VDJ', 'VH'])
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
