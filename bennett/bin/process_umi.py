#!/usr/bin/env python
# encoding: utf-8
"""
"""
import sys
import gzip
import editdist as ed
import subprocess as sp
from toolshed import nopen
from collections import Counter, defaultdict
from itertools import islice, groupby, izip

__version__ = "0.2"

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
    cmd += ("""awk '{printf($0);n++;if(n%4==0){printf("\\n")}else{printf("\\t")}}' | """
                """awk '{i=substr($2,1,length); print i"\\t"$0}' | """
                """sort -k1,1 -k2,2 | """
                """cut -f 2,3,4,5 | """
                """tr "\\t" "\\n\"""").replace("length", args.length)
    sp.call(cmd, shell=True)

def decode(x):
    """
    >>> decode("4")
    19
    """
    return ord(x) - 33

def average(quals):
    """
    >>> average("/==96996<FGCHHHGGGFFE=EDFFFEEB")
    32.53...
    """
    vals = map(decode, quals)
    return sum(vals)/float(len(quals))

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


def process_pairs(r1out, r2out, r1seqs, r2seqs, r1seq_to_name, r2name_to_seq, r1seq_to_qual, r2seq_to_qual, umi, read_id):
    """prints most abundant sequence observed across a particular UMI.
    
    r1out - file handle
    r2out - file handle
    r1seqs - counter for R1
    r2seqs - counter for R2
    r1seq_to_name - dictionary to lookup name by sequence
    r2name_to_seq - dictionary to lookup sequence by name
    r1seq_to_qual - dictionary to obtain quality score
    r2seq_to_qual - dictionary to obtain quality score
    umi - the sequence being processed
    read_id - digit for unique output read ids
    """
    # most_common_r1_seq, most_common_r1_count = r1seqs.most_common(1)[0]
    # most_common_r2_seq, most_common_r2_count = r2seqs.most_common(1)[0]
    # 
    # print >>sys.stderr, most_common_r1_count
    # print >>sys.stderr, most_common_r2_count
    score_cutoff = 4
    best_score = 0
    r1_score = 0
    r2_score = 0
    best_r1_seq = ""
    best_r2_seq = ""

    # i think 3 is enough tries; leaving it hardcoded
    for r1_sequence, r1_count in r1seqs.most_common(3):
        r2seq_counter = Counter()
        for name in r1seq_to_name[r1_sequence]:
            r2seq_counter.update([r2name_to_seq[name]])
        r2_most_common_seq, r2_most_common_count = r2seq_counter.most_common(1)[0]
        scores = [r1_count, r2_most_common_count]
        score = sum(scores) * (min(scores) / float(max(scores)))
        
        # print >>sys.stderr, "{score}: {r1} {r2}".format(score=score, r1=r1_count, r2=r2_most_common_count)
        
        if score > best_score:
            best_score = score
            best_r1_seq = r1_sequence
            r1_score = r1_count
            best_r2_seq = r2_most_common_seq
            r2_score = r2_most_common_count

    if best_score >= score_cutoff:
        print >>sys.stderr, "{best}: {r1} {r2}\n".format(best=best_score, r1=r1_score, r2=r2_score)
        r1out.write("@r_%d:%s 1\n%s\n+\n%s\n" % (read_id, umi, best_r1_seq, r1seq_to_qual[best_r1_seq]))
        r2out.write("@r_%d:%s 2\n%s\n+\n%s\n" % (read_id, umi, best_r2_seq, r2seq_to_qual[best_r2_seq]))
        read_id += 1
    return read_id

def get_name(name, insert):
    sample = name.split(".")[0]
    return "{sample}.{insert}.fastq.gz".format(sample=sample, insert=insert)

def add_qual(d, k, q):
    try:
        if average(d[k]) < average(q):
            d[k] = q
    except KeyError:
        d[k] = q
    return d

def run_collapse(args):
    """r1, r2, umi, mismatches"""
    mmatch = args.mismatches
    iupac_umi = args.umi
    umileng = len(iupac_umi)
    cutoff = args.cutoff
    readid = 1
    r1out = gzip.open(get_name(args.r1, "umifiltered"), 'wb')
    r2out = gzip.open(get_name(args.r2, "umifiltered"), 'wb')

    for umi, group in groupby(izip(readfq(args.r1), readfq(args.r2)), key=lambda (rr1, rr2): rr1.seq[:umileng]):
        if not valid_umi(iupac_umi, umi): continue
        r1seqs = Counter()
        r1seq_to_name = defaultdict(list)
        r2name_to_seq = {}
        r1seq_to_qual = {}
        r2seq_to_qual = {}

        r2seqs = Counter()

        for r1, r2 in group:
            assert r2.seq[:umileng] == umi
            assert r1.name.split()[0] == r2.name.split()[0]
            
            # trim UMI and clip at first N
            trimmed_r1_seq = r1.seq.split("N", 1)[0][umileng:]
            trimmed_r2_seq = r2.seq.split("N", 1)[0][umileng:]

            if len(trimmed_r1_seq) < cutoff or len(trimmed_r2_seq) < cutoff: continue

            trimmed_r1_qual = r1.qual[umileng:len(trimmed_r1_seq) + umileng]
            trimmed_r2_qual = r2.qual[umileng:len(trimmed_r2_seq) + umileng]
            
            assert len(trimmed_r1_seq) == len(trimmed_r1_qual)
            assert len(trimmed_r2_seq) == len(trimmed_r2_qual)

            r1seqs.update([trimmed_r1_seq])
            r1seq_to_name[trimmed_r1_seq].append(r1.name)
            r2name_to_seq[r2.name] = trimmed_r2_seq
            
            r2seqs.update([trimmed_r2_seq])
            
            # maintains best qual per seq
            r1seq_to_qual = add_qual(r1seq_to_qual, trimmed_r1_seq, trimmed_r1_qual)
            r2seq_to_qual = add_qual(r2seq_to_qual, trimmed_r2_seq, trimmed_r2_qual)

        # process UMI group, writing to files, and returning current read id
        readid = process_pairs(r1out, r2out, r1seqs, r2seqs, r1seq_to_name, r2name_to_seq, r1seq_to_qual, r2seq_to_qual, umi, readid)

def main(args):
    args.func(args)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, version=__version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subp = p.add_subparsers(help='commands')

    # sorting the fastq by umi
    fsort = subp.add_parser('sort', description="Sorts fastq by UMI.", help="order the fastq by the UMI to facilitate processing")
    fsort.add_argument('fastq', metavar='FASTQ', help='unsorted reads with incorporated UMI')
    fsort.add_argument('length', metavar='LENGTH', help='length of the UMI sequence')
    fsort.set_defaults(func=run_sort)
    
    # pull out each unique sequence per UMI
    fcollapse = subp.add_parser('collapse', description="Finds unique sequences per valid UMI among paired-end reads.", help="find most abundant sequence per UMI given paired-end reads")
    fcollapse.add_argument('r1', metavar="R1", help="R1 FASTQ with UMI, sorted by UMI.")
    fcollapse.add_argument('r2', metavar="R2", help="R2 FASTQ with UMI, sorted by UMI.")
    fcollapse.add_argument('umi', metavar="UMI", help='IUPAC sequence of the UMI, e.g. NNNNNV')
    fcollapse.add_argument('-c', '--cutoff', type=int, default=200, help='shortest allowable read length after trimming at first N')
    fcollapse.add_argument('-m', '--mismatches', type=int, default=3, help='allowable mismatches when finding unique sequences')
    fcollapse.set_defaults(func=run_collapse)
    
    args = p.parse_args()
    main(args)
