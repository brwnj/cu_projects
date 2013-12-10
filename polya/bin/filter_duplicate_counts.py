#!/usr/bin/env python
# encoding: utf-8
"""
filter the site counts removing entries where gene is accounted for with same
count and same peak classification, but has multiple refseq annotations.

NM_181886	1.NM_181886|ube2d3.2	3
NM_181886	1.NM_181886|ube2d3.4	2
NM_181886	1.NM_181886|ube2d3.5	3
NM_181886	1.NM_181886|ube2d3.6	3
NM_181886	1.NM_181886|ube2d3.7	2
NM_181887	1.NM_181887|ube2d3.2	3
NM_181887	1.NM_181887|ube2d3.4	2
NM_181887	1.NM_181887|ube2d3.5	3
NM_181887	1.NM_181887|ube2d3.6	3
NM_181887	1.NM_181887|ube2d3.7	2
NM_181888	1.NM_181888|ube2d3.2	3
NM_181888	1.NM_181888|ube2d3.4	2
NM_181888	1.NM_181888|ube2d3.5	3
NM_181888	1.NM_181888|ube2d3.6	3
NM_181888	1.NM_181888|ube2d3.7	2
NM_181889	1.NM_181889|ube2d3.2	3
NM_181889	1.NM_181889|ube2d3.4	2
NM_181889	1.NM_181889|ube2d3.5	3
NM_181889	1.NM_181889|ube2d3.6	3
NM_181889	1.NM_181889|ube2d3.7	2
NM_181890	1.NM_181890|ube2d3.2	3
NM_181890	1.NM_181890|ube2d3.4	2
NM_181890	1.NM_181890|ube2d3.5	3
NM_181890	1.NM_181890|ube2d3.6	3
NM_181890	1.NM_181890|ube2d3.7	2
NM_181891	1.NM_181891|ube2d3.2	3
NM_181891	1.NM_181891|ube2d3.4	2
NM_181891	1.NM_181891|ube2d3.5	3
NM_181891	1.NM_181891|ube2d3.6	3
NM_181891	1.NM_181891|ube2d3.7	2
NM_181892	1.NM_181892|ube2d3.2	3
NM_181892	1.NM_181892|ube2d3.4	2
NM_181892	1.NM_181892|ube2d3.5	3
NM_181892	1.NM_181892|ube2d3.6	3
NM_181892	1.NM_181892|ube2d3.7	2
NM_181893	1.NM_181893|ube2d3.2	3
NM_181893	1.NM_181893|ube2d3.4	2
NM_181893	1.NM_181893|ube2d3.5	3
NM_181893	1.NM_181893|ube2d3.6	3
NM_181893	1.NM_181893|ube2d3.7	2
"""
from toolshed import reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def main(count_file):
    observed = {}
    for toks in reader(count_file, header=['gene', 'fullname', 'count']):
        if not "|" in toks['fullname']: continue
        # ube2d3.2
        gene_symbol = toks['fullname'].split("|")[1]
        peak_class = toks['fullname'].split(".", 1)[0]


        if observed.has_key(gene_symbol):
            try:
                if observed[gene_symbol][peak_class] = toks['count']: continue
            except KeyError:
                observed[gene_symbol]


        print "\t".join(t for t in toks.values())




if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    p.add_argument(counts)
    args = p.parse_args()
    main(args.counts)
