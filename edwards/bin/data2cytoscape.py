#!/usr/bin/env python
# encoding: utf-8
"""
parse data table into cytoscape files
"""

from toolshed import reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def main(table):
    d = {}
    for toks in reader(table, header=True, sep=" "):
        row_gene = toks['Genes']
        d[row_gene] = {}
        for col_gene in toks.keys():
            # row 1, col 1 is a generic header entry
            if col_gene == "Genes": continue
            d[row_gene][col_gene] = int(toks[col_gene])

    # print node size attributes
    node_out = open("node_attrs.txt", "wb")
    print >>node_out, "source\ttotal_mutations"
    for k in d.keys():
        print >>node_out, "{gene}\t{count}".format(gene=k, count=d[k][k])
    node_out.close()

    # print network and edge attributes
    interaction_type = "pp"
    network_out = open("network.txt", "wb")
    print >>network_out, "source\tinteraction_type\ttarget\tcomutation_count"
    seen = set()
    for row_gene in d.keys():
        for col_gene, count in d[row_gene].iteritems():
            if count == 0: continue
            # double checking these were filtered out
            if row_gene == col_gene: continue
            # check to see if the interaction was already added in the opposite direction
            if "{gene2}_{gene1}".format(gene2=col_gene, gene1=row_gene) in seen: continue
            print >>network_out, "{gene1}\t{interaction}\t{gene2}\t{count}".format(gene1=row_gene, interaction=interaction_type, gene2=col_gene, count=count)
            seen.add("{gene1}_{gene2}".format(gene1=row_gene, gene2=col_gene))
    network_out.close()

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('table')
    args = p.parse_args()
    main(args.table)
