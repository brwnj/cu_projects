#!/usr/bin/env python
# encoding: utf-8
"""
Add gene descriptions onto snpeff output.
"""
import argparse
from toolshed import reader, header

def descriptions_to_dict(fname):
    d = dict()
    for l in reader(fname):
        d[l['systematic_name']] = l
    return d

def main(snpeff, descriptions):
    snpeff_header = header(snpeff)
    additional_columns = header(descriptions)
    full_header = snpeff_header + additional_columns
    descriptions_dict = descriptions_to_dict(descriptions)
    print "\t".join(full_header)
    for l in reader(snpeff):
        #Gene_ID Gene_name
        # l['Gene_ID']
        add = dict()
        try:
            add = descriptions_dict[l['Gene_ID']]
        except KeyError:
            try:
                add = descriptions_dict[l['Gene_name']]
            except KeyError:
                print "\t".join(l[h] for h in snpeff_header)
                continue
        full_line = dict(l.items() + add.items())
        try:
            assert full_line.get('status')
        except AssertionError:
            full_line['status'] = ""
        print "\t".join(full_line[h] for h in full_header)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("snpeff_output")
    p.add_argument("gene_descriptions")
    args = p.parse_args()
    main(args.snpeff_output, args.gene_descriptions)
