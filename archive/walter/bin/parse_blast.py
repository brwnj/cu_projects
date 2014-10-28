#!/usr/bin/env python
"""
Parse blastn output. Find alignments that map in opposite orientations.

Example output:
chr1:1000273-1000313	ENST00000334019_+	93.75	16	1	0	12	27	211	196	  127	24.3

From command:
blastn -query peaks.fa -db ensg_3utr -evalue 1000 -word_size 7 -gapopen 5 -gapextend 2 -reward 1 -penalty -3 -outfmt 6

-outfmt 6 by default returns:
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

"""
import sys
from toolshed import reader

def main(args):
    blast_output = "qseqid sseqid pident length mismatch gapopen qstart qend \
                    sstart send evalue bitscore".split()
    # might need to add more stuff
    add_output = "genecard".split()
    # print a file header
    blast_output.extend(add_output)
    print "\t".join(blast_output)
    for b in reader(args.blast, header=blast_output):
        # parse qseqid for proper strand orientation [+, -]
        ensembl_id, known_strand = b['sseqid'].split("_")
        b['genecard'] = "some url for gene cards"
        # determine strand alignment
        aligned_strand = "+" if int(b['sstart']) < int(b['send']) else "-"
        if known_strand != aligned_strand:
            print "\t".join([b[i] for i in blast_output])

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('blast', help='std tabular blast output')
    main(p.parse_args())