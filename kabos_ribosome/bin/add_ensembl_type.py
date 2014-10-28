#!/usr/bin/env python
# coding=utf-8
"""
Append type category onto DESeq output from Ensembl's gtf.
"""

import os
import pickle

from pybedtools import BedTool
from toolshed import reader, header
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def make_ens_table(gtf):
    t = {}
    gtf = BedTool(gtf)
    for toks in gtf:
        t[toks.attrs['gene_id']] = toks[1]
    return t


def main(deseq_output, ensembl_gtf):

    p = "%s.pickle" % ensembl_gtf
    if os.path.exists(p):
        with open(p, 'rb') as handle:
            ensembl_translation = pickle.load(handle)
    else:
        ensembl_translation = make_ens_table(ensembl_gtf)
        with open(p, 'wb') as handle:
            pickle.dump(ensembl_translation, handle)

    h = header(deseq_output, sep=",")
    h[0] = 'ensg'
    h = [i.strip('"') for i in h]
    h.append('type')

    print ",".join(h)

    for toks in reader(deseq_output, sep=","):
        toks['ensg'] = toks['']
        toks['type'] = ensembl_translation[toks['ensembl']]
        print ",".join(toks[i] for i in h)


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('deseq_output', help="csv with ensg as first column")
    p.add_argument('ensembl_gtf', help="self-explanatory")
    args = vars(p.parse_args())
    main(**args)
