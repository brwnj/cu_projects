#!/usr/bin/env python
# encoding: utf-8
"""
Add gene descriptions onto snpeff output.
"""
import os
import tempfile
from toolshed import header, nopen, reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def fix_header(snpeff):
    f = open(tempfile.mkstemp(suffix=".txt")[1], "w")
    for line in nopen(snpeff):
        line = line.strip()
        if line.startswith("#"):
            toks = line.split()
            if toks[1] == "Chromo":
                line = line.lstrip("# ")
                print >>f, line
        else:
            print >>f, line
    f.close()
    return f.name

def descriptions_to_dict(fname):
    d = {}
    for l in reader(fname):
        d[l['Ensembl Gene ID']] = l
    return d

def main(snpeff, descriptions):
    fixed_header = fix_header(snpeff)
    snpeff_header = header(fixed_header)
    additional_columns = header(descriptions)
    full_header = snpeff_header + additional_columns
    descriptions_dict = descriptions_to_dict(descriptions)
    print "\t".join(full_header)
    for l in reader(fixed_header):

        # lines from snpeff have varying number of fields
        # need to fill in blanks
        for tok in snpeff_header:
            if l.has_key(tok): continue
            l[tok] = "na"

        add = {}

        try:
            add = descriptions_dict[l['Gene_ID']]
        except KeyError:
            try:
                add = descriptions_dict[l['Gene_name']]
            except KeyError:
                print "\t".join(l[h] for h in snpeff_header)
                continue

        full_line = dict(l.items() + add.items())
        print "\t".join(full_line[h] for h in full_header)

    os.remove(fixed_header)

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("snpeff_output")
    p.add_argument("gene_descriptions")
    args = p.parse_args()
    main(args.snpeff_output, args.gene_descriptions)
