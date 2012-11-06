#!/usr/bin/env python
# encoding: utf-8
"""
Writes a count matrix to stdout. Requires the raw counts of individual samples
and the ENSEMBL reference gene model.

File name should include the ID like: xx_ID.xx

Created by Joe Brown on 2012-02-15.
"""
import argparse
import os
import sys
import re
from toolshed import nopen


def makeref(gtf):
    """takes a gtf and creates a dictionary -- ENSG:GENENAME"""
    sys.stderr.write("Preparing gene name cross-reference.\n")
    ref = {}
    previousgene = ''
    for line in nopen(gtf):
        ens = re.findall(r'gene_id \"([\w+\.]+)\"', line)[0].strip()
        if ens == previousgene:
            continue
        previousgene = ens
        name = re.findall(r'gene_name \"([a-zA-Z-_\/|.0-9]+)\"', line)[0].strip()
        ref[ens] = name
    return ref


def countsmatrix(gtf, files):
    """Prints count data to stdout."""
    generef = makeref(gtf)
    sys.stderr.write("Preparing counts table.\n")
    all_data = {}
    for f in files:
        sample = os.path.basename(f).rsplit(".")[0].split("_")[1]
        all_data[sample] = {}
        for line in open(f):
            toks = line.split("\t")
            #htseq annotated lines of filtered genes
            if not toks[0].startswith("ENSG"):
                continue
            name = '%s,%s' % (toks[0], generef.get(toks[0]))
            count = str(int(toks[1]))
            all_data[sample][name] = count
    samples = sorted(all_data)
    genes = all_data[samples[1]].keys()
    print "#ID\t" + "\t".join(samples)
    for gene in genes:
        print gene + "\t" + "\t".join(all_data[s][gene] for s in samples)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    required = p.add_argument_group('required arguments')
    required.add_argument("-gtf", dest="gtf", help="ENSEMBL GTF", type=str)
    required.add_argument("-f",dest="files",nargs="+",help="Count files from htseq-count.")
    args = p.parse_args()
    if args.gtf and args.files:
        countsmatrix(args.gtf, args.files)
    else:
        p.print_help()
        sys.exit(1)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
