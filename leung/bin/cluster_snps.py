#!/usr/bin/env python
# coding=utf-8
"""
"""
import sys
import pandas as pd
from toolshed import reader
from collections import Counter
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def main(files):
    c = {}
    for txt in files:
        sample = txt.split("_annotated")[0]
        c[sample] = Counter()
        for toks in reader(txt, header=True, skip_while=lambda x: not x[0].startswith('# Chromo')):
            if not toks['Effect'] == "NON_SYNONYMOUS_CODING"\
                or not int(toks['Coverage']) > 10\
                or not toks['Change_type'] == "SNP": continue

            c[sample].update([toks['Gene_ID']])

    df = pd.DataFrame(c)
    df.to_csv("clustered_genes.csv", na_rep=0)
    df.T.to_csv("clustered_genes_transposed.csv", na_rep=0)


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('files')
    args = vars(p.parse_args())
    main(**args)
