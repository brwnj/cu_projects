#!/usr/bin/env python
# encoding: utf-8
"""
Finding unique restriction site associated DNA (RAD) tag sequences per unique
molecular identifier (UMI).
"""
import sys
import editdist as ed
import subprocess as sp
from toolshed import nopen
from collections import Counter, defaultdict
from itertools import islice, groupby, izip

__version__ = "1.0"

IUPAC = {"A":"A","T":"T","C":"C","G":"G","R":"GA","Y":"TC",
         "M":"AC","K":"GT","S":"GC","W":"AT","H":"ACT",
         "B":"GTC","V":"GCA","D":"GAT","N":"GATC"}

class Fastq(object):
    def __init__(self, args):
        self.name = args[0][1:]
        self.seq = args[1]
        self.qual = args[3]
        assert len(self.seq) == len(self.qual)
    
    def __repr__(self):
        return "Fastq({name})".format(name=self.name)
    
    def __str__(self):
        return "@{name}\n{seq}\n+\n{qual}".format(name=self.name,
                seq=self.seq, qual=self.qual)

def readfq(fq):
    with nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            rd = [x for x in islice(fqclean, 4)]
            if not rd: raise StopIteration
            assert all(rd) and len(rd) == 4
            yield Fastq(rd)

def valid_umi(iupac, umi):
    """parse UMI sequence to validate against IUPAC sequence."""
    for code, base in izip(iupac, umi):
        try:
            if not base in IUPAC[code]:
                return False
        except KeyError:
            return False
    return True

def run_sort(args):
    """fastq, length"""
    if args.fastq.endswith(".gz"):
        cmd = "gunzip -c %s | " % args.fastq
    else:
        cmd = "cat %s | " % args.fastq
    # sort by the UMI, then by read name for paired-end sync
    cmd += """awk '{printf($0);n++;if(n%4==0){printf("\\n")}else{printf("\\t")}}' |\
                awk '{i=substr($2,1,8); print i"\\t"$0}' | sort -k1,1 -k2,2 |\
                cut -f 2,3,4,5 | tr "\\t" "\\n\""""
    sp.call(cmd, shell=True)

def trimmed_seq(seq, leng, t5, t3):
    """docstring for trimmed_seq"""
    return seq[leng + t5:-t3] if t3 > 0 else seq[leng + t5:]

def distance(a, b):
    """find best edit distance between two strings of potentially uneven length.
    """
    la, lb = len(a), len(b)
    if la < lb:
        return distance(b, a)
    if la == lb:
        return ed.distance(a, b)
    else:
        dists = []
        for i in xrange(0, la-lb+1):
            dists.append(ed.distance(a[i:i+lb], b))
        return min(dists)

def process_pairs(r1_out, r2_out, r1counter, r1dd, r2d, umiseq, mismatches, readid):
    """prints every unique read per umi"""
    ignore = set()
    seen = set()
    r1seqs = []
    r1seqs = set(list(r1counter))
    for target in r1seqs:
        if target in seen: continue
        seen.add(target)
        for query in r1seqs:
            if query in seen: continue
            if distance(target, query) < mismatches:
                for name in r1dd[query]:
                    # add similar read names to appropriate bin
                    r1dd[target].append(name)
                ignore.add(query)
                seen.add(query)
    chosen_seqs = r1seqs - ignore
    # find most abundant R2 and print reads
    for seq in chosen_seqs:
        r2seqs = Counter()
        for name in r1dd[seq]:
            # list of read names used in this bin
            r2seqs.update([r2d[name]])
        # print the fasta records
        r1_out.write(">read_%d:%s 1\n%s\n" % (readid, umiseq, seq))
        r2_out.write(">read_%d:%s 2\n%s\n" % (readid, umiseq, \
                                                r2seqs.most_common(1)[0][0]))
        readid += 1
    return readid

def run_scanp(args):
    """r1i, r2i, r1o, r2o, umi, mismatches, r1five, r1three, r2five, r2three"""
    leng = len(args.umi)
    readid = 0
    with nopen(args.r1i) as r1i, nopen(args.r2i) as r2i,\
            open(args.r1o, 'w') as r1o, open(args.r2o, 'w') as r2o:
        r1iter = read_fastq(r1i)
        r2iter = read_fastq(r2i)
        for (r1name, r1seq, r1qual), (r2name, r2seq, r2qual) in izip(r1iter, r2iter):

            umiseq = r1seq[:leng]
            r1name = r1name.split()[0]
            r2name = r2name.split()[0]
            # really no longer need the UMI on R2
            assert umiseq == r2seq[:leng]
            assert r1name == r2name
            # skip anything with an invalid UMI
            if not valid_umi(args.umi, umiseq): continue

            # count unique sequences
            r1seqs = Counter()
            # sequence to read name
            r1dd = defaultdict(list)
            # add the current sequence
            r1seq_trim = trimmed_seq(r1seq, leng, args.r1five, args.r1three)
            r2seq_trim = trimmed_seq(r2seq, leng, args.r2five, args.r2three)

            r1seqs.update([r1seq_trim])
            r1dd[r1seq_trim].append(r1name)
            # dictionary of name to sequence
            r2d = {r2name:r2seq_trim}

            # get the remaining sequences of this UMI
            for match, group in groupby(izip(r1iter, r2iter), \
                    lambda ((n1,s1,q1),(n2,s2,q2)): s1[:leng] == umiseq):
                if match:
                    # compile all the reads of this UMI using
                    for (g1name, g1seq, g1qual), (g2name, g2seq, g2qual) in group:
                        g1name = g1name.split()[0]
                        g2name = g2name.split()[0]
                        g1seq_trim = trimmed_seq(g1seq, leng, \
                                                    args.r1five, args.r1three)
                        g2seq_trim = trimmed_seq(g2seq, leng, \
                                                    args.r2five, args.r2three)
                        # counter
                        r1seqs.update([g1seq_trim])
                        # defaultdict(list)
                        r1dd[g1seq_trim].append(g1name)
                        # dictionary
                        r2d[g2name] = g2seq_trim
                else:
                    break
            readid = process_pairs(r1o, r2o, r1seqs, r1dd, r2d, \
                                    umiseq, args.mismatches, readid)

def main(args):
    args.func(args)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, version=__version__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    subp = p.add_subparsers(help='commands')

    # sorting the fastq by umi
    fsort = subp.add_parser('sort',
            description="Sorts fastq by UMI.",
            help="order the fastq by the UMI to facilitate processing")
    fsort.add_argument('fastq', metavar='FASTQ',
            help='unsorted reads with incorporated UMI')
    fsort.add_argument('length', metavar='LENGTH',
            help='length of the UMI sequence')
    fsort.set_defaults(func=run_sort)
    
    # find most abundance sequence per umi among paired-end reads
    fscanp = subp.add_parser('scanp',
            description="Finds most abundant sequence per valid UMI among \
                    paired-end reads.",
            help="find most abundant sequence per UMI given paired-end reads")
    fscanp.add_argument('r1i', metavar="R1-IN",
            help="R1 FASTQ with UMI to scan.")
    fscanp.add_argument('r2i', metavar="R2-IN",
            help="R2 FASTQ with UMI to scan.")
    fscanp.add_argument('r1o', metavar="R1-OUT", help="R1 output FASTQ.")
    fscanp.add_argument('r2o', metavar="R2-OUT", help="R2 output FASTQ.")
    fscanp.add_argument('umi', metavar="UMI",
            help='IUPAC sequence of the UMI, e.g. NNNNNV')
    fscanp.add_argument('-m', '--mismatches', type=int, default=3,
            help='allowable mismatches when finding unique sequences [%(default)s]')
    fscanp.set_defaults(func=run_scanp)
    
    args = p.parse_args()
    main(args)
