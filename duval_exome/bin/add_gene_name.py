#!/usr/bin/env python
# encoding: utf-8
"""
Add gene_name to canFam3 GTF from UCSC.
"""

import re
from toolshed import nopen, reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def main(gtf, xref):
    xref_dict = {}
    for toks in reader(xref, header=['ens', 'symbol']):
        xref_dict[toks['ens']] = toks['symbol']

    gene_id = re.compile(r"gene_id \"(.*?)\"")
    for line in nopen(gtf):
        line = line.rstrip("\r\n")
        toks = line.split("\t")
        annotations = toks[-1]
        ens = gene_id.findall(toks[-1])[0]
        assert ens.startswith("ENS")
        try:
            print line + ' gene_name "' + xref_dict[ens] + '";'
        except KeyError:
            print line + ' gene_name "' + ens + '";'

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__,
                formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('gtf')
    p.add_argument('xref')
    args = p.parse_args()
    main(args.gtf, args.xref)
