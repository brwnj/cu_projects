#!/usr/bin/env python
# encoding: utf-8
"""
Add gene descriptions onto snpeff output.
"""
from toolshed import reader, header
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def descriptions_to_dict(fname):
    d = {}
    for l in reader(fname):
        d[l['Ensembl Gene ID']] = l
    return d

def main(snpeff, descriptions):
    snpeff_header = header(snpeff)
    additional_columns = header(descriptions)
    full_header = snpeff_header + additional_columns
    descriptions_dict = descriptions_to_dict(descriptions)
    print "\t".join(full_header)
    for l in reader(snpeff):
        if not l.has_key(snpeff_header[0]): continue
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

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("snpeff_output")
    p.add_argument("gene_descriptions")
    args = p.parse_args()
    main(args.snpeff_output, args.gene_descriptions)
