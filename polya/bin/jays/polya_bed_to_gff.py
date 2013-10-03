#! /usr/bin/env python

import sys
from tabdelim import DictReader

# GTF2 format from http://mblab.wustl.edu/GTF2.html
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

# example BED line
# chr1  227718952   227718953   p.c1.ABCB10.3   0   -

source = 'polya'
fieldnames = ('chrom','chromStart','chromEnd','name','score','strand')
bedreader = DictReader(file(sys.argv[1]), fieldnames=fieldnames)

for row in bedreader:
    seqname = row['name']
    start = int(row['chromStart']) + 1 # 1-based in GTF
    end = int(row['chromEnd']) + 1 # 1-based in GTF
    score = row['name'].split('.')[1].replace('c','')
    feature = 'pA.class.%s' % score
    strand = row['strand']
    frame = 0
    group = row['name'].split('.')[2]

    fields = (seqname, source, feature, start, end,
              score, strand, frame, group)
    print '\t'.join(map(str, fields))

