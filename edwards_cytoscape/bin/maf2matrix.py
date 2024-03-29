#!/usr/bin/env python
# coding=utf-8
"""
Parsing .maf
"""

import pandas as pd
from collections import Counter
from itertools import combinations
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def main(data_file, gene, effect, vtype, patient, find_effect, find_type, network_file, min_comutation, nodeattrs_file):
    df = pd.read_table(data_file, compression="gzip" if data_file.endswith(".gz") else None)
    # filter out irrelevant mutations
    df = df[(df[effect].isin(find_effect)) & (df[vtype].isin(find_type))]

    interactions = Counter()
    for _, gdf in df.groupby(patient):
        # all genes for an individual
        genes = gdf[gene].values.tolist()
        seen = set()
        # all pair combinations for genes list
        for (a, b) in combinations(genes, 2):
            # add pairs to counter
            # but only count a gene:gene interaction once per patient
            pair = "%s:%s" % (a, b)
            if pair in seen: continue
            seen.add(pair)
            interactions.update([pair])

    with open(network_file, 'w') as out:
        interaction_type = "pp"
        # add header onto file
        print >>out, "\t".join(['source', 'interaction_type', 'target', 'comutation_count'])
        for gene_pair, comutation_count in interactions.iteritems():
            if comutation_count < min_comutation: continue
            source, target = gene_pair.split(":")
            if source == target: continue
            print >>out, "%s\t%s\t%s\t%d" % (source, interaction_type, target, comutation_count)

    total_counts = Counter()
    # this is clearly not the best way of doing this...
    for (g, s), count in df.groupby([gene, patient]).size().iteritems():
        total_counts.update([g])
    with open(nodeattrs_file, 'w') as out:
        print >>out, "\t".join(['source', 'total_mutations'])
        for gene_name, count in total_counts.iteritems():
            print >>out, "%s\t%d" % (gene_name, count)


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('data_file', help="maf data file with header")

    data_file_args = p.add_argument_group("required metadata")
    data_file_args.add_argument('--gene', default="Hugo_Symbol", help="column name containing gene ID")
    data_file_args.add_argument('--effect', default="Variant_Classification", help="column name containing variant effect")
    data_file_args.add_argument('--vtype', default="Variant_Type", help="column name containing variant type")
    data_file_args.add_argument('--patient', default="Tumor_Sample_Barcode", help="column name containing sample identification")

    user_specs_args = p.add_argument_group("parse parameters")
    user_specs_args.add_argument('-e', '--find-effect', default=["Missense_Mutation"], action="append", help="mutation effect to quantify; can be specified multiple times")
    user_specs_args.add_argument('-t', '--find-type', default=["SNP"], action="append", help="mutation type to quantify; can be specified multiple times")

    output_file_args = p.add_argument_group("output files")
    output_file_args.add_argument('--network-file', default="network.txt", help="file name containing network file")
    output_file_args.add_argument('--min-comutation', default=2, type=int, help="minimum comutation count in order to output interaction")
    output_file_args.add_argument('--nodeattrs-file', default="node_attrs.txt", help="file name containing node attributes")

    args = vars(p.parse_args())
    main(**args)
