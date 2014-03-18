#!/usr/bin/env python
# encoding: utf-8
"""
Filter reads mapped from trimmed fastq where the tag was not identified in the
read, but is located within `bases` of `discard`.
"""

from os import remove
from sys import stderr, exit
from pysam import Samfile, index
from subprocess import check_call
from os.path import splitext, exists
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def bed_from_bam(bam, mark, n):
    """
    returns name of bed file.

    bam  - has aligned reads marked with `mark` that we need regions for
    mark - unique string in the read name
    n    - number of bases downstream to expand the region
    """
    # index the bam file if not found
    if not exists(bam + ".bai"):
        print >>stderr, "indexing", bam
        index(bam)

    total_marked = 0
    bed = open(splitext(bam)[0] + "_bed_regions.temp", 'w')
    with Samfile(bam) as bam_file:
        for chrom in bam_file.references:
            for read in bam_file.fetch(chrom):
                if read.is_unmapped: continue

                if mark == 'all' or mark in read.qname:
                    total_marked += 1
                    # coordinates will include read sequence in fasta
                    if read.is_reverse:
                        start = read.pos - n if read.pos > n else 1
                        end = read.pos + read.alen
                        strand = "-"
                    else:
                        start = read.pos
                        end = start + read.alen + n
                        strand = "+"

                    # name has to include position since we're allowing reads
                    # to map in more than one location
                    b = [chrom, start, end, read.qname + ":" + str(read.pos), "0", strand]
                    print >>bed, "\t".join([str(f) for f in b])
    bed.close()

    if total_marked == 0:
        print >>stderr, "nothing was found to be removed with:", mark
        remove(bed.name)
        exit(1)

    return bed.name, total_marked


def seqs_from_bed(bed, fasta):
    """
    returns name of seqs file.

    bed   - intervals to pull sequences
    fasta - reference fasta for species
    """
    seqs = bed.split("_bed_regions.temp")[0] + "_seqs.temp"
    try:
        cmd = ["bedtools", "getfasta", "-s", "-tab", "-name", "-fi",
                fasta, "-bed", bed, "-fo", seqs]
        check_call(cmd)
    finally:
        remove(bed)
    return seqs


def find_discards(sequences, discard):
    """
    returns set of read names.

    sequences - files of read name and sequence with tab in between
    discard   - string of bases; if found read will be removed from bam
    """
    discards = set()
    try:
        with open(sequences) as fh:
            for line in fh:
                name, sequence = line.strip().split("\t")

                for dseq in discard:

                    if sequence.find(dseq) == -1: continue
                    discards.add(name)
                    break

    finally:
        remove(sequences)
    return discards


def filter_bam(inbam, outbam, discards):
    """
    returns number of filtered reads.

    inbam    - unfiltered mapped reads
    outbam   - filtered mapped reads
    discards - set of read ids to remove from inbam
    """
    filtered_reads = 0
    with Samfile(inbam, 'rb') as in_bam, \
            Samfile(outbam, 'wb', template=in_bam) as out_bam:
        for read in in_bam.fetch():
            if read.qname + ":" + str(read.pos) in discards:
                filtered_reads += 1
            else:
                out_bam.write(read)

    return filtered_reads


def main(inbam, outbam, fasta, bases, discard, mark):
    assert exists(inbam)
    assert exists(fasta)

    print >>stderr, ">> finding marked reads"
    bed, total_marked = bed_from_bam(inbam, mark, bases)
    print >>stderr, ">> pulling sequences for", total_marked, "locations"
    sequences = seqs_from_bed(bed, fasta)
    print >>stderr, ">> analyzing downstream sequences"
    discards = find_discards(sequences, discard)
    print >>stderr, ">> writing filtered bam"
    filtered_reads = filter_bam(inbam, outbam, discards)
    # very slight discordance between filtered_reads and len(discards)
    # MCF10A -- 384311, 384222; not going to investigate now
    print >>stderr, "Filtered", filtered_reads, "of", total_marked, "marked for inspection"


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('inbam', help="bam to filter")
    p.add_argument('outbam', help="filtered bam")
    p.add_argument('fasta', help="reference fasta; index required")
    p.add_argument('-b', '--bases', default=100, type=int, help="number of bases downstream to look for discard sequence")
    p.add_argument('-d', '--discard', default=["TCAGTC"], action="append", help="if found without flanking adapter piece, discard mapped read")
    p.add_argument('-m', '--mark', default="inspect", help="string that was used to mark reads or 'all' if you want to check all aligned reads")

    args = vars(p.parse_args())
    main(**args)
