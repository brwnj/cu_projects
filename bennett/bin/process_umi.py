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

__version__ = "0.1"

IUPAC = {"A":"A","T":"T","C":"C","G":"G","R":"GA","Y":"TC",
         "M":"AC","K":"GT","S":"GC","W":"AT","H":"ACT",
         "B":"GTC","V":"GCA","D":"GAT","N":"GATC"}

def readfq(fq):
    class Fastq(object):
        def __init__(self, args):
            self.name = args[0][1:].partition(" ")[0]
            self.seq = args[1]
            self.qual = args[3]
            assert len(self.seq) == len(self.qual)
    
        def __repr__(self):
            return "Fastq({name})".format(name=self.name)
    
        def __str__(self):
            return "@{name}\n{seq}\n+\n{qual}".format(name=self.name,
                    seq=self.seq, qual=self.qual)
    
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

# def trimmed_seq(seq, leng, t5, t3):
#     return seq[leng + t5:-t3] if t3 > 0 else seq[leng + t5:]

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

def process_pairs(r1out, r2out, r1seqs, r1seq_to_name, r2name_to_seq, umi, m, n):
    """prints every unique read per umi and add the UMI sequence to the name."""
    ignore = set()
    seen = set()
    r1seqsset = set(list(r1seqs))
    for target in r1seqsset:
        if target in seen: continue
        seen.add(target)
        for query in r1seqsset:
            if query in seen: continue
            if distance(target, query) < m:
                for name in r1seq_to_name[query]:
                    # add similar read names to appropriate bin
                    r1seq_to_name[target].append(name)
                ignore.add(query)
                seen.add(query)
    chosen_seqs = r1seqsset - ignore
    # find most abundant R2 and print reads
    for seq in chosen_seqs:
        r2seqs = Counter()
        for name in r1seq_to_name[seq]:
            # list of read names used in this bin
            r2seqs.update([r2name_to_seq[name]])
        
        # need a fastq
        r1out.write(">read_%d:%s 1\n%s\n" % (readid, umiseq, seq))
        r2out.write(">read_%d:%s 2\n%s\n" % (readid, umiseq, r2seqs.most_common(1)[0][0]))
        n += 1
    return n

def get_name(name):
    if ".fastq" in name:
        sample = name.split(".fastq")[0]
    else:
        sample = name.split(".fq")[0]
    return "{sample}.umifiltered.fastq".format(**locals())

def run_scanp(args):
    """r1, r2, umi, mismatches"""
    mmatch = args.mismatches
    iupac_umi = args.umi
    leng = len(iupac_umi)
    readid = 1
    r1out = open(get_name(args.r1), 'wb')
    r2out = open(get_name(args.r2), 'wb')
    for umi, group in groupby(izip(readfq(args.r1), readfq(args.r2)), key=lambda (rr1, rr2): rr1.seq[:leng]):
        r1seqs = Counter()
        r1seq_to_name = defaultdict(list)
        r2name_to_seq = {}
        for r1, r2 in group:
            assert r2.seq[:leng] == umi
            assert r1.name.split()[0] == r2.name.split()[0]
            if not valid_umi(iupac_umi, umi): continue
            
            trimmed_r2_seq = r2.seq.split("N")[0]
            if len(trimmed_r2_seq) < 100: continue
            trimmed_r2_qual = r2.qual[len(trimmed_r2_seq):]

            r1seqs.update([r1.seq])
            r1seq_to_name[r1.seq].append(r1.name)
            r2name_to_seq[r2.name] = trimmed_r2_seq
        
        # magic...
        readid = process_pairs(r1out, r2out, r1seqs, r1seq_to_name, r2name_to_seq, umi, mmatch, readid)
        sys.exit(1)

def main(args):
    args.func(args)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, version=__version__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subp = p.add_subparsers(help='commands')

    # sorting the fastq by umi
    fsort = subp.add_parser('sort', description="Sorts fastq by UMI.", help="order the fastq by the UMI to facilitate processing")
    fsort.add_argument('fastq', metavar='FASTQ', help='unsorted reads with incorporated UMI')
    fsort.add_argument('length', metavar='LENGTH', help='length of the UMI sequence')
    fsort.set_defaults(func=run_sort)
    
    # pull out each unique sequence per UMI
    fscanp = subp.add_parser('scanp', description="Finds unique sequences per valid UMI among paired-end reads.", help="find most abundant sequence per UMI given paired-end reads")
    fscanp.add_argument('r1', metavar="R1", help="R1 FASTQ with UMI to scan.")
    fscanp.add_argument('r2', metavar="R2", help="R2 FASTQ with UMI to scan.")
    fscanp.add_argument('umi', metavar="UMI", help='IUPAC sequence of the UMI, e.g. NNNNNV')
    fscanp.add_argument('-m', '--mismatches', type=int, default=3, help='allowable mismatches when finding unique sequences [%(default)s]')
    fscanp.set_defaults(func=run_scanp)
    
    args = p.parse_args()
    main(args)
