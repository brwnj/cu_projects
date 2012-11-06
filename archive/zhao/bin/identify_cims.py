#! /usr/bin/env python

'''
cliptools-identify-cims: identify CIMS (deletions) within collapsed reads.

N.B. because all reads from a given alignment need to be analyzed to produce an
accurate number, this should not be parallelized.

alignments must have been run such that deletions are reported in the
alignment report. (crossmatch "-discrep_tags" option).
'''

import sys
import pysam
from collections import defaultdict, Counter
from functools import partial

from seqtools.formats import read_crossmatch
from genomedata._util import maybe_gzip_open

__version__ = '$Revision: 404 $'
__author__ = 'Jay Hesselberth <jay.hesselberth@gmail.com>'

# TODO:
#  - add support for novoalign

def identify_cims(align_filenames, verbose):

    # counter of deletions observed at specific sites across
    # unique reads
    deletions = defaultdict(Counter)

    # keep track of which reads have deletions in them so we only count
    # them once
    defaultset = partial(defaultdict, set)
    seen_deletions = defaultdict(defaultset)

    for filename in align_filenames:

        if verbose:
            msg = ">> identifying CIMS in: %s"
            print msg % filename
            
        if filename.endswith(".bam"):
            with pysam.Samfile(filename, "rb") as align_file:
                for alignedread in align_file:
                    #chrom =
                    pass
        elif filename.endswith(".sam"):
            with pysam.Samfile(filename, "r") as align_file:
                for alignedread in align_file:
                    pass
        else:                                                                   #don't like this at all since this is calling per read discrepancies
            with maybe_gzip_open(filename) as align_file:
                for record in read_crossmatch(align_file):

                    chrom = record.alignment.target_id
                    start = record.alignment.target_start
                    end = record.alignment.target_end
                    strand = crossmatch_strand(record)

                    # parse the descrepancies
                    for discrep in record.discreps:

                        # XXX need to check whether stranded data is being
                        # handled correctly
                        # TODO improve this by adding option to write out data
                        # by mutation type (D,I or S)

                        # need to include target_pos and discrep.type
                        # in info as some reads will have adjacent mutations 
                        target_pos = int(discrep.target_pos)
                        site_info = (end, strand, target_pos, discrep.type)

                        # check whether an overlapping read has already been
                        # observed to have a deletion. must do the check here
                        # as there could be other overlapping reads that do
                        # NOT have deletions, and thus would never be counted
                        # if the check is done above.
                        if not site_info in seen_deletions[chrom][start]:
                            deletions[chrom][target_pos] += 1

                        # mark the deletion-containing read as seen
                        seen_deletions[chrom][start].add(site_info)

    # write out the deletions in BEDGRAPH format
    for chrom in sorted(deletions):
        for start, count in sorted(deletions[chrom].items()):
            fields = (chrom, start, start + 1, count)
            print '\t'.join(map(str, fields))

def crossmatch_strand(record):
    if record.complement:
        return '-'
    return '+'

def parse_options(args):
    from optparse import OptionParser

    description = ("Identify CIMS sites in alignment files")
    usage = '%prog [options] ALIGNFILES...'
    version = '%%prog %s' % __version__

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    parser.add_option("-v", "--verbose",
        action="store_true", dest="verbose",
        help="maximum verbosity [default: %default]",
        default=False)

    options, args = parser.parse_args(args)

    if len(args) == 0:
        parser.error("specify alignment files")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    align_filenames = args
    kwargs = {'verbose':options.verbose}

    return identify_cims(align_filenames, **kwargs)

if __name__ == '__main__':
    sys.exit(main())