#!/usr/bin/env python
# encoding: utf-8
"""
Returns an unsorted bed from UCSC.

Dependencies:
pymysql

Created by Joe Brown on 2011-12-21.
"""
import argparse
import os
import sys
try:
    import pymysql
except ImportError:
    sys.stderr.write("pymsql not found.\nTry: pip install pymsql\n")
    sys.exit(1)


def getbed(user, host, species, port, table, strip):
    """Returns a bed from UCSC tables."""
    conn = pymysql.connect(host=host, port=port, user=user, db=species)
    cur = conn.cursor()
    if table is 'ensGene':
        cur.execute("select chrom, txStart, txEnd, T.name, \
                    X.value as geneName, strand, exonStarts, \
                    exonEnds, S.source from %s as T, ensemblToGeneName as X, \
                    ensemblSource as S where X.name=T.name and S.name=T.name"
                    % table)
    else:
        cur.execute("select chrom, txStart, txEnd, T.name, \
                    X.geneSymbol as geneName, strand, exonStarts, \
                    exonEnds from %s as T, kgXref as X \
                    where X.kgID=T.name" %
                    table)
    header = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                'exonstart', 'exonend', 'color', 'count', 'sizes', 'starts']
    for r in cur.fetchall():
        line = {}
        if strip:
            line['chrom'] = r[0]
        else:
            line['chrom'] = r[0].lstrip("chr")
        line['start'] = str(r[1])
        line['end'] = str(r[2])
        if table is 'ensGene':
            line['name'] = '%s;%s;%s' % (r[3], r[4], r[8]) # extra biotype col
        else:
            line['name'] = '%s;%s' % (r[3], r[4])
        line['score'] = '1'
        line['strand'] = r[5]
        line['exonstart'] = str(r[1])
        line['exonend'] = str(r[2])
        line['color'] = '.'
        astarts = r[6].split(',')
        aends = r[7].split(',')
        starts = []
        sizes = []
        count = 0
        for i in range(len(astarts)):
            if not astarts[i]: continue
            sizes.append('%d' % (int(aends[i]) - int(astarts[i])))
            starts.append('%d' % (int(astarts[i]) - r[1]))
            count = count + 1
        line['count'] = str(count)
        line['sizes'] = ','.join(s for s in sizes)
        line['starts'] = ','.join(s for s in starts)
        print "\t".join(line[h] for h in header)
    cur.close()
    conn.close()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    required = p.add_argument_group('required arguments')
    p.add_argument("--no-chr-strip", dest='strip', action='store_true', 
        default=False, help="Set to retain 'chr' on chromosome name.")
    p.add_argument("--species", dest='species', nargs='?', default='hg19')
    p.add_argument("--table", dest='table', nargs='?', default='ensGene', 
        help="Other possible option is 'knownGene'")
    p.add_argument("--user", dest="user", nargs='?', default='genome')
    p.add_argument("--host", dest="host", nargs='?', 
        default='genome-mysql.cse.ucsc.edu')
    p.add_argument("--port", dest="port", nargs='?', default='3306', type=int)
    args = p.parse_args()
    getbed(args.user, args.host, args.species, args.port, args.table, args.strip)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()