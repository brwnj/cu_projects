#!/usr/bin/env python
# encoding: utf-8
"""
Join R1 and R2, or print unpaired reads of R1 or R2

This script assumes there is whitespace separating the readname from the pair number:

                                       > <
@HW-ST997:140:D0UA5ACXX:8:1101:1497:2138 1:N:0:ATCACG
@HW-ST997:140:D0UA5ACXX:8:1101:1497:2138 2:N:0:ATCACG
                                       > <
"""
import sys
import bz2
import gzip
import os.path as op

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def fastqparser(fastq):
    """yields name, seq, qual from fastq file"""
    line_num = -1
    record = []
    for line in fopen(fastq):
        line_num += 1
        if line_num == 4:
            yield record[0][1:], record[1], record[3]
            line_num = 0
            record = []
        record.append(line.strip())    
    if record:
        if record[0]:
            yield record[0][1:], record[1], record[3]


def fopen(f, mode="rb"):
    """open file, use gzip.open if gzipped."""
    if f.endswith((".gz", ".Z", ".z")):
        return gzip.open(f, mode)
    elif f.endswith((".bz", ".bz2", ".bzip2")):
        return bz2.BZ2File(f, mode)
    else:
        return open(op.expanduser(op.expandvars(f)), mode)


def reverse(s): 
    """return in reverse order""" 
    seq = list(s) 
    seq.reverse() 
    return ''.join(seq)


def complement(s):
    """return complementary seq""" 
    basecomplement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
                        'a':'t', 'c':'g', 'g':'c', 't':'a'} 
    seq = list(s)
    seq = [basecomplement[b] for b in seq]
    return ''.join(seq)


def reversecomplement(s):
    """return reverse complement"""
    s = reverse(s)
    s = complement(s)
    return s


def fastqtodict(fastq):
    """returns dict of read name to sequence"""
    fdict = {}
    for name, seq, qual in fastqparser(fastq):
        # explicitly state space to facilitate future changes
        fdict[name.split(" ")[0]] = seq
    return fdict


def main(args):
    # set dictionary based on mode
    if args.mode == "R2":
        fastqdict = fastqtodict(args.R1)
        fastq = args.R2
    else:
        fastqdict = fastqtodict(args.R2)
        fastq = args.R1
    
    for name, seq, qual in fastqparser(fastq):
        try:
            # explicitly state space to facilitate future changes
            name = name.split(" ")[0]
            cseq = fastqdict.get(name)
            if args.reverse_complement:
                cseq = reversecomplement(cseq)
            print "@%s\n%s%s" % (name, seq, cseq)
        except KeyError:
            # without pairs
            if not args.mode == "paired":
                print "@%s\n%s" % (name, seq)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('R1', help="fastq of R1")
    p.add_argument('R2', help="fastq of R2")
    p.add_argument("--mode", "-m", choices=["paired", "R1", "R2"], 
                    default="paired", help="prints <mode> to stdout (default: paired)")
    p.add_argument("--reverse_complement", "-r", action="store_false", 
                    default=True, help="reverse complement R2 when joining (default: True)")
    args = p.parse_args()
    main(args)