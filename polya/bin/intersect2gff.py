#!/usr/bin/env python
# encoding: utf-8
"""
bedtools intersect output to dexseq_counts.py gff format
"""
import re
# import sys
from toolshed import reader
from collections import defaultdict, OrderedDict

def main(args):
    genes = defaultdict(OrderedDict)
    # the s is for "shit i don't need"
    for toks in reader(args.intersect, header="chrom s s start stop s strand s attrs s s s genename s s".split()):
        polyaname = re.findall(r'transcript_id \"([a-zA-Z-_\/|.0-9]+)\"', toks['attrs'])[0].strip()
        transcript = {}
        transcript[polyaname] = {'chrom':toks['chrom'],'start':toks['start'],'stop':toks['stop'],'strand':toks['strand']}
        genes[toks['genename']].update(transcript)

    for gene, transcripts in genes.iteritems():
        source = "intersect2gff"
        genestart = genestop = None
        chrom = ""
        strand = ""
        # get start and stop of the gene
        for polya, toks in transcripts.iteritems():
            chrom = toks['chrom']
            strand = toks['strand']
            if genestart is None:
                genestart = int(toks['start'])
            if genestart > int(toks['start']):
                genestart = int(toks['start'])
            if genestop is None:
                genestop = int(toks['stop'])
            if genestop < int(toks['stop']):
                genestop = int(toks['stop'])
        # print gene line
        genefields = [chrom, source, "aggregate_gene", genestart - 1, genestop + 1, ".", strand, ".", 'gene_id "%s"' % gene]
        print "\t".join(map(str, genefields))
        # print each polya site
        # counter = 0
        for polya, toks in transcripts.iteritems():
            # counter += 1
            # attrs = 'transcripts "%s"; exonic_part_number "%03d"; gene_id "%s"' % (polya, counter, gene)
            attrs = 'transcripts "%s"; exonic_part_number "%s"; gene_id "%s"' % (polya, polya, gene)
            polyafields = [toks['chrom'], source, "exonic_part", int(toks['start']) - 1, int(toks['stop']) + 1, ".", toks['strand'], ".", attrs]
            print "\t".join(map(str, polyafields))

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('intersect', help='result file of `bedtools intersect -wb -s -a polya_3utronly.gtf -b refseq.utr3.bed`')
    args = p.parse_args()
    main(args)