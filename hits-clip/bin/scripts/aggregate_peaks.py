#!/usr/bin/env python

'''cliptools-aggregate-peaks: aggregate peaks by presence of seed
sequences and report parent genes. '''

import sys
import string
from toolshed import reader, nopen
from genomedata import Genome
from collections import defaultdict, Counter

__version__ = '1.0'
__author__ = 'Jay Hesselberth; brwnj'


def revcomp(s, _comp=string.maketrans('ATCG', 'TAGC')):
    return s.translate(_comp)[::-1]


def aggregate_peaks(peak_bedfilename, gene_bedfilename, seed_fastafilename,
                    gdfilename, strand, spec_chrom, verbose):

    seeds = load_seeds(seed_fastafilename, verbose)
    peaks = load_peaks(peak_bedfilename, verbose)
    genes = load_genes(gene_bedfilename, verbose)

    if not strand:
        strand = find_strand_from_filename(peak_bedfilename)
        if verbose:
            print >>sys.stderr, ">> learned strand %s from filename" % \
                strand

    with Genome(gdfilename) as genome:

        if verbose:
            print >>sys.stderr, ">> analyzing peaks..."

        for peak_chrom in peaks:

            if spec_chrom and peak_chrom != spec_chrom: continue

            for peak_start in peaks[peak_chrom]:

                peak_end = peaks[peak_chrom][peak_start]

                chrom = genome[peak_chrom]

                # want to compare the reverse complement of the peak
                # sequence to the seed, where the seqs should be the same
                try:
                    peak_seq = chrom.seq[peak_start:peak_end].tostring()
                except NotImplementedError:
                    msg = ">> no seq available for %s:%d-%d"
                    print >>sys.stderr, msg % (peak_chrom, peak_start,
                                               peak_end)
                    continue

                if strand == '+':
                    peak_seq = revcomp(peak_seq)

                for seed_id, seed_seq in seeds.items():
                    if seed_seq in str(peak_seq):
                        parent_gene = find_parent_gene(genes, peak_chrom, peak_start,
                                                       peak_end)

                        if not parent_gene: continue

                        fields = (seed_id, peak_chrom, peak_start, peak_end,
                                  parent_gene)
                        print '\t'.join([str(i) for i in fields])


def find_strand_from_filename(filename):
    if 'pos' in filename:
        return '+'
    elif 'neg' in filename:
        return '-'
    else:
        raise TypeError, "Cannot learn strand from bed filename"


def find_parent_gene(genes, chrom, start, end):
    for gene_start in genes[chrom]:
        if start >= gene_start:
            gene_end, gene_name, gene_score, gene_strand = \
                genes[chrom][gene_start]
            if end <= gene_end:
                return gene_name

    return None


def readfa(fa):
    with nopen(fa) as fh:
        for header, group in groupby(fh, lambda line: line[0] == '>'):
            if header:
                line = group.next()
                name = line[1:].strip()
            else:
                seq = ''.join(line.strip() for line in group)
                yield name, seq


def load_seeds(seed_fastafilename, verbose):
    if verbose:
        print >>sys.stderr, ">> loading seeds..."

    seeds = {}
    for name, seq in readfa(seed_fastafile):
        seeds[name] = seq

    return seeds


def load_peaks(peak_bedfilename, verbose):

    if verbose:
        print >>sys.stderr, ">> loading peaks..."

    peaks = defaultdict(dict)

    for peak in reader(peak_bedfilename, header="chrom start stop name score strand".split()):
            data = int(peak['stop'])
            peaks[peak['chrom'][int(peak['start']] = data

    return peaks


def load_genes(gene_bedfilename, verbose):

    if verbose:
        print >>sys.stderr, ">> loading genes..."

    genes = defaultdict(dict)

    for gene in reader(gene_bedfilename, header="chrom start stop name score strand".split()):
        data = (int(gene['stop']), gene['name'], int(gene['score']), gene['strand'])
        genes[gene['chrom']][int(gene['start'])] = data

    return genes


def parse_options(args):
    from optparse import OptionParser

    description = ("Identify miRNA seeds in peak regions.")
    usage = '%prog [options] PEAK_BED GENE_BED SEED_FASTA GDFILENAME'
    version = '%%prog %s' % __version__

    parser = OptionParser(usage=usage, version=version,
                          description=description)


    parser.add_option("-c", "--chrom",
        action="store", dest="chrom",
        help="peak chrom (default: %default)",
        default=None)

    parser.add_option("-s", "--strand",
        action="store", dest="strand",
        help="peak strand (default: %default)",
        default=None)

    parser.add_option("-v", "--verbose",
        action="store_true", dest="verbose",
        help="maximum verbosity (default: %default)",
        default=False)

    options, args = parser.parse_args(args)

    if not len(args) == 4:
        parser.error("specify PEAKS, GENES, SEEDS and GENOMEDATADIR")

    return options, args

def main(args=sys.argv[1:]):
    options, args = parse_options(args)

    peak_bedfilename, gene_bedfilename, seed_fastafilename, \
        gdfilename = args
    kwargs = {'verbose':options.verbose,
              'strand':options.strand,
              'spec_chrom':options.chrom}

    return aggregate_peaks(peak_bedfilename, gene_bedfilename,
        seed_fastafilename, gdfilename, **kwargs)

if __name__ == '__main__':
    sys.exit(main())

